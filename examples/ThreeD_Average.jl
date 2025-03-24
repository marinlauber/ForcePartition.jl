using WaterLily,StaticArrays,WriteVTK,JLD2
using ForcePartition

function circ(D,dims,Ubc;Re=200,U=1,mem=CUDA.CuArray,perdir=())
    body = AutoBody((x,t)->√sum(abs2,SA[x[1],x[2]].-dims[2]÷2)-D÷2)
    Simulation(dims, Ubc, D; body, ν=U*D/Re, mem, perdir)
end

import WaterLily: scale_u!,BCTuple,conv_diff!,accelerate!,BDIM!,BC!,project!,CFL,@loop,size_u
@fastmath function mom_step_force!(a::Flow{N},b::AbstractPoisson,∇τ) where N
    a.u⁰ .= a.u; scale_u!(a,0)
    # predictor u → u'
    U = BCTuple(a.U,@view(a.Δt[1:end-1]),N)
    conv_diff!(a.f,a.u⁰,a.σ,ν=a.ν,perdir=a.perdir)
    accelerate!(a.f,@view(a.Δt[1:end-1]),a.g,a.U)
    force!(a,∇τ) # new forcing
    BDIM!(a); BC!(a.u,U,a.exitBC,a.perdir)
    a.exitBC && exitBC!(a.u,a.u⁰,U,a.Δt[end]) # convective exit
    project!(a,b); BC!(a.u,U,a.exitBC,a.perdir)
    # corrector u → u¹
    U = BCTuple(a.U,a.Δt,N)
    conv_diff!(a.f,a.u,a.σ,ν=a.ν,perdir=a.perdir)
    accelerate!(a.f,a.Δt,a.g,a.U)
    force!(a,∇τ) # new forcing
    BDIM!(a); scale_u!(a,0.5); BC!(a.u,U,a.exitBC,a.perdir)
    project!(a,b,0.5); BC!(a.u,U,a.exitBC,a.perdir)
    push!(a.Δt,CFL(a))
end
@inline force!(a::Flow,Ⅎ::AbstractArray) = @loop a.f[Ii] += Ⅎ[Ii] over Ii in CartesianIndices(a.f)

# make a writer with some attributes, need to output to CPU array to save file (|> Array)
_velocity(a::Simulation) = a.flow.u |> Array;
_pressure(a::Simulation) = a.flow.p |> Array;
vort(a::Simulation) = (@loop a.flow.f[I,:] .= WaterLily.ω(I,a.flow.u) over I ∈ inside(a.flow.p);
                       a.flow.f |> Array)
_body(a::Simulation) = (measure_sdf!(a.flow.σ, a.body, WaterLily.time(a)); 
                                     a.flow.σ |> Array;)
lamda(a::Simulation) = (@loop a.flow.σ[I] = WaterLily.λ₂(I, a.flow.u) over I ∈ inside(a.flow.p);
                        a.flow.σ |> Array;)

function u_span(a::Simulation) # spanwise average velocity
    ϵ = inv(size(inside(a.flow.p),3)) # sum over domain
    N,n = size_u(a.flow.u); a.flow.f .= 0 # reset
    for i ∈ 1:2
        @loop a.flow.f[Base.front(I),1,i] = ϵ*a.flow.u[I,i] over I ∈ inside(a.flow.p)
        @loop a.flow.f[I,i] = a.flow.f[Base.front(I),1,i] over I ∈ inside(a.flow.p)
    end
    return a.flow.f |> Array
end

custom_attrib = Dict(
    "u" => _velocity,
    "p" => _pressure,
    "ω" => vort,
    "b" => _body,
    "λ₂" => lamda,
    "us" => u_span
)# this maps what to write to the name in the file
# make the writer
using CUDA
R = 16 # 3*2^5
# star with 2D sim, run it to get wake, then run 3D and SANS version
sim2D = circ(2R,(20R,8R),(1,0);Re=10_000,U=1,mem=CuArray)
sim_step!(sim2D,50;verbose=true) # 2D update

# full 3D, initialize with 2D sim
sim = circ(2R,(20R,8R,2R),(1,0,0);Re=10_000,U=1,mem=CuArray,perdir=(3,))
spread!(sim2D.flow,sim.flow; ϵ=1e-3) # spread 2D to 3D
avrg = SpanAverage(sim.flow) # span average structure

# writers
writer = vtkWriter("ThreeD_SANS"; attrib=custom_attrib)
writer2D = vtkWriter("ThreeD_SANS_2D")

# reset 2D sim time
foreach(i->pop!(sim2D.flow.Δt),1:length(sim2D.flow.Δt)-1)

# data = jldopen("WaterLily_3D_SANS.jld2")
# sim.flow.u .= data["u"]
# sim.flow.p .= data["p"]
# span_average!(avrg,sim.flow)
# SANS!(avrg.v,avrg.τ)

for t in range(0,20.0;step=0.05)#1:6
    while sim_time(sim)<t #sim_step!(sim,t)
        mom_step!(sim.flow,sim.pois) # 3D update
        span_average!(avrg,sim.flow)
        SANS!(avrg.v,avrg.τ) # compute forcing ∇⋅τ
        mom_step_force!(sim2D.flow,sim2D.pois,avrg.v) # 2D update
        sim2D.flow.Δt[end] = sim.flow.Δt[end] # keep time in sync
    end
    write!(writer,sim); write!(writer2D,sim2D)
    @show t
    flush(stdout)
end
# jldsave("WaterLily_3D_SANS.jld2";
#     u=Array(sim.flow.u),
#     p=Array(sim.flow.p),
#     t=time(sim)
# )
close(writer); close(writer2D)