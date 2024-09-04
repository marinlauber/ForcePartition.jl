using WaterLily
using Plots
using StaticArrays
using ForcePartition
using ParametricBodies
using BiotSavartBCs
include("TwoD_plots.jl")

function make_sphere(;L=32,Re=250,U=1)
    radius, center = L/2, 2L
    center = SA[center,center,center]

    # make a body
    sphere = AutoBody((x,t)->√sum(abs2, x .- center .- 1.5) - radius)

    # generate sim
    Simulation((8L,4L,4L), (U,0,0), radius; ν=U*radius/Re, body=sphere)
end

# make simulation
sim = make_sphere(L=32)

# contruct FPM and compute potential
fpm = ForcePartitionMethod(sim)
potential!(fpm;tᵢ=0,f=:force,axis=1,itmx=32)

# evolve the simulation
sim_step!(sim,1)
@inside sim.flow.σ[I] = WaterLily.curl(3,I,sim.flow.u)*sim.L/sim.U
flood(sim.flow.σ[:,:,2sim.L],clims=(-5,5))

# qfield influence
qϕ = -∫2Qϕ!(fpm,sim.flow,20.0,recompute=false)
# total force
force = -WaterLily.pressure_force(sim)

# run a sim and plot the time evolution
sim = make_sphere(L=32)
fpm = ForcePartitionMethod(sim)
t₀,duration,step = 0.,10,0.2
force, fp = [],[]

@time @gif for tᵢ in range(t₀,t₀+duration;step)
    # update until time tᵢ in the background
    while sim_time(sim) < tᵢ
        # update flow
        mom_step!(sim.flow,sim.pois)
        # pressure force
        push!(force,-2WaterLily.pressure_force(sim)[1])
        push!(fp,-∫2Qϕ!(fpm,sim.flow,recompute=false))
    end

    # plot vorticity
    @inside sim.flow.σ[I] = WaterLily.curl(3,I,sim.flow.u)*sim.L/sim.U
    flood(sim.flow.σ; shift=(-0.5,-0.5),clims=(-5,5))
    body_plot!(sim); plot!(title="tU/L $tᵢ")

    # print time step
    println("tU/L=",round(tᵢ,digits=4),", Δt=",round(sim.flow.Δt[end],digits=3))
end

# # plot results
# plot(cumsum(@views(sim.flow.Δt[1:end-1]))./sim.L,force,label="-2∫pndS",xlabel="tU/L",ylabel="Force")
# plot!(cumsum(@views(sim.flow.Δt[1:end-1]))./sim.L,fp,label="-2∫QϕdV")
