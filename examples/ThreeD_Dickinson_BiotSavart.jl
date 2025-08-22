using WaterLily,StaticArrays,BiotSavartBCs
using CUDA,ForcePartition,WriteVTK#,Plots

# make a writer with some attributes, need to output to CPU array to save file (|> Array)
vtk_velocity(a::AbstractSimulation) = a.flow.u |> Array;
vtk_pressure(a::AbstractSimulation) = a.flow.p |> Array;
vtk_body(a::AbstractSimulation) = (measure_sdf!(a.flow.σ, a.body, WaterLily.time(a)); a.flow.σ |> Array;)
vtk_lambda(a::AbstractSimulation) = (@inside a.flow.σ[I] = WaterLily.λ₂(I, a.flow.u); a.flow.σ |> Array;)
vtk_Q(a::AbstractSimulation) = (@inside a.flow.σ[I] = ForcePartition.Qcriterion(I, a.flow.u); a.flow.σ |> Array)
vtk_ϕ(a::AbstractSimulation) = (potential!(fpm,a.body;tᵢ=sum(a.flow.Δt)); fpm.ϕ |> Array)

# strating from rest
WaterLily.CFL(a::Flow) = WaterLily.CFL(a;Δt_max=0.5)

# overwrite to apply symmetric BiotBCs
import BiotSavartBCs: interaction,_interaction!
@inline function _interaction!(ml,lT)
    (l,T) = lT     # level & target
    n = length(ml) # total levels
    ω = ml[l]      # vorticity sources

    T⁺,sgn = image(T,size(ω),2) # image target & sign
    ml[l][T] = interaction(ω,T,l,n)+sgn*interaction(ω,T⁺,l,n)
end

# make simulation following Dickinson's setup
function Dickinson(L=64;U=1,Re=5e2,ε=0.5f0,thk=2ε+√3,AR=2,mem=Array,T=Float32)

    # ellipse sdf
    function sdf(x,t)
        # offset it to put origin at tip, and stretch it by AR. remove radius for sdf
        √sum(abs2, SA[x[1], (x[2]-1.5f0L/2.f0)/AR, 0]) - L/2/AR
    end

    # the mapping
    function map(x,t)
        # Dickinson kinemtics
        _α = π/2.f0 - π/4.f0*sin(π*t/L)  # positive pitch increase AoA
        _ϕ = 0.35f0π*cos(π*t/L)      # positive to the rear of the mosquito
        # rotation mmatrix
        Ry = SA[cos(_α) 0 sin(_α); 0 1 0; -sin(_α) 0 cos(_α)] # alpha
        Rz = SA[cos(_ϕ) -sin(_ϕ) 0; sin(_ϕ) cos(_ϕ) 0; 0 0 1] # phi
        return Ry*Rz*(x .- SA[2L,0,L]) # the order matters
    end

    # Build the mosquito from a mapped ellipsoid and two plane that trim it to the correct thickness
    ellipsoid = AutoBody(sdf, map)
    upper_lower = AutoBody((x,t)->(abs(x[3])-thk/2.f0), map)
    body = ellipsoid ∩ upper_lower # intersection of sets

    # Return initialized simulation
    return BiotSimulation((4L,3L,4L),(0,0,0),L;ν=U*L/Re,U,body,mem,T,nonbiotfaces=(2,))
end

# make a simulation
sim = Dickinson(96;mem=CuArray,T=Float32)

# make a vtk writer
custom_attrib = Dict("u" => vtk_velocity,"p" => vtk_pressure,"d" => vtk_body,
                     "λ2" => vtk_lambda) #,"Q" => vtk_Q,"ϕ" => vtk_ϕ)# this maps what to write to the name in the file
# make the writer
writer = vtkWriter("Dickinson"; attrib=custom_attrib)

# force moment
# fpm = ForcePartitionMethod(sim,axis=3) # axis to 3 for lift
# potential!(fpm,fpm.body)

# a simulation
t₀ = sim_time(sim); duration = 12; tstep = 0.05 # print time
forces,fp,fm,fv = [],[],[],[] # empty list to store forces

# simulate
@time for tᵢ in range(t₀,t₀+duration;step=tstep)
    # update until time tᵢ in the background
    t = sum(sim.flow.Δt[1:end-1])
    # this is normaly hiden into sim_step!(sim,tᵢ)
    while t < tᵢ*sim.L/sim.U
        #update the body
        sim_step!(sim;remasure=true)
        # update time
        t += sim.flow.Δt[end]
        # compute and save pressure forces, this is actually a force coefficient,
        # to find the actual force we have Fi = Ci/(1/2*ρ*U^2*L^2)
        force = -2WaterLily.pressure_force(sim)
        push!(forces,[t/sim.L,force...])
        # the totential is the same, we just recompute it once and then use it
        # push!(fp,-∫2QϕdV!(fpm,sim.flow,recompute=true))
        # push!(fm,-∮UϕdS!(fpm,sim.flow,recompute=false))
        # push!(fv, ∮ReωdS!(fpm,sim.flow,recompute=false))
    end
    # this writes every tstep to paraview files, not every time step
    save!(writer,sim);
    println("tU/L=",round(tᵢ,digits=4),", Δt=",round(sim.flow.Δt[end],digits=3))
end
forces = reduce(vcat,forces')
println("Done...")
close(writer)

# # plot results
# time = cumsum(sim.flow.Δt)./sim.L # time vector
# plot(time[1:end-1],[forces[:,2]./sim.L^2, fp/sim.L^2, fm/sim.L^2, fv/sim.L^2, (fp+fm+fv)/sim.L^2],
#      label=["Total force" "vorticity force" "added-mass force" "viscous force" "total partition"],
#      xlabel="tU/L",ylabel="2F/ρU²L²") #,ylims=(1,3),xlims=(0,12))
# savefig("force_partition_mosquito.png")