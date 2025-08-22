using WaterLily,StaticArrays,WriteVTK,JLD2

function circ(D,dims,Ubc;Re=200,U=1,mem=CUDA.CuArray,perdir=())
    body = AutoBody((x,t)->√sum(abs2,SA[x[1],x[2]].-dims[2]÷2)-D÷2)
    Simulation(dims, Ubc, D; body, ν=U*D/Re, mem, perdir)
end

import WaterLily: @loop,inside_u,size_u
# f = ∇⋅τ = ∂ᵢτᵢⱼ, for example f₁ = ∂₁τ₁₁ + ∂₂τ₂₁, because τ is at cell center inline components are
# simply the difference between neighboring cells, but the diagonal components are interpolated using the
# average of the two neighboring cells.
@fastmath SANS!(f::AbstractArray, τ::AbstractArray) = (f.=0; for (i,j) ∈ Iterators.product(1:2,1:2)
    @loop f[I,i] = (i==j ? τ[I,i,j]-τ[I-δ(i,I),i,j] :
                    @inbounds(τ[I+δ(j,I),i,j]+τ[I-δ(i,I)+δ(j,I),i,j]
                             -τ[I-δ(j,I),i,j]-τ[I-δ(i,I)-δ(j,I),i,j])/4)  over I in inside_u(f)
end)

# add the SANS terms
@inline force!(a::Flow, t; Ⅎ) = @loop a.f[Ii] += Ⅎ[Ii] over Ii in CartesianIndices(a.f)

# make a writer with some attributes, need to output to CPU array to save file (|> Array)
vtk_velocity(a::AbstractSimulation) = a.flow.u |> Array;
vtk_pressure(a::AbstractSimulation) = a.flow.p |> Array;
vtk_vort(a::AbstractSimulation) = (@loop a.flow.f[I,:] .= WaterLily.ω(I,a.flow.u) over I ∈ inside(a.flow.p);
                                    a.flow.f |> Array)
vtk_body(a::AbstractSimulation) = (measure_sdf!(a.flow.σ, a.body, WaterLily.time(a));
                                    a.flow.σ |> Array;)
vtk_lamda(a::AbstractSimulation) = (@loop a.flow.σ[I] = WaterLily.λ₂(I, a.flow.u) over I ∈ inside(a.flow.p);
                                    a.flow.σ |> Array;)

function vtk_u_span(a::AbstractSimulation) # spanwise average velocity
    ϵ = inv(size(inside(a.flow.p),3)) # sum over domain
    N,n = size_u(a.flow.u); a.flow.f .= 0 # reset
    for i ∈ 1:2
        @loop a.flow.f[Base.front(I),1,i] = ϵ*a.flow.u[I,i] over I ∈ inside(a.flow.p)
        @loop a.flow.f[I,i] = a.flow.f[Base.front(I),1,i] over I ∈ inside(a.flow.p)
    end
    return a.flow.f |> Array
end

custom_attrib = Dict("u"=>vtk_velocity, "p"=>vtk_pressure, "ω"=>vtk_vort,
                     "b" => vtk_body, "λ₂" => vtk_lamda, "us"=>vtk_u_span)# this maps what to write to the name in the file

# make the writer
using CUDA
R = 16 # 3*2^5
# star with 2D sim, run it to get wake, then run 3D and SANS version
sim2D = circ(2R,(20R,8R),(1,0);Re=10_000,U=1,mem=CuArray)
sim_step!(sim2D,5;verbose=true) # 2D update

# full 3D, initialize with 2D sim
sim = circ(2R,(20R,8R,2R),(1,0,0);Re=10_000,U=1,mem=CuArray,perdir=(3,))
WaterLily.spread!(sim2D.flow,sim.flow; ϵ=1e-3) # spread 2D to 3D
avrg = SpanAverage(sim.flow) # span average structure
update!(avrg, sim.flow) # update average with initial condition

# writers
writer = vtkWriter("ThreeD_SANS"; attrib=custom_attrib)
writer2D = vtkWriter("ThreeD_SANS_2D")

# reset 2D sim time
foreach(i->pop!(sim2D.flow.Δt), 1:length(sim2D.flow.Δt)-1)

Ⅎ = avrg.U; # pointer

for t in range(sim_time(sim),sim_time(sim)+20.0;step=0.05)
    while sim_time(sim)<t
        sim_step!(sim) # 3D update
        update!(avrg, sim.flow)
        SANS!(Ⅎ, avrg.UU) # compute forcing ∇⋅τ
        sim_step!(sim2D; udf=force!, Ⅎ) # 2D update
        sim2D.flow.Δt[end] = sim.flow.Δt[end] # keep time in sync
    end
    save!(writer,sim); save!(writer2D,sim2D)
    @show t
    flush(stdout)
end
# jldsave("WaterLily_3D_SANS.jld2";
#     u=Array(sim.flow.u),
#     p=Array(sim.flow.p),
#     t=time(sim)
# )
close(writer); close(writer2D)