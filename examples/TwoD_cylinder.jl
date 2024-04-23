using WaterLily
using Plots
using StaticArrays
using ForcePartition
include("TwoD_plots.jl")

function make_sim(L;Re=250,U=1)
    radius, center = L/8, L/2
    center = SA[center,center]

    # make a body
    circle = AutoBody((x,t)->√sum(abs2, x .- center .- 1.5) - radius)

    # generate sim
    Simulation((2*L,L), (U,0), radius; ν=U*radius/Re, body=circle)
end

# make simulation
sim = make_sim(2^8)

# contruct fpm and compute potential
fpm = ForcePartitionMethod(sim)
potential!(fpm,0)

# plot to check
R = inside(fpm.ϕ)
flood(fpm.ϕ[R],clims=(-32,32),levels=11,filled=true)
body_plot!(sim)

# evolve the simulation
sim_step!(sim,10;verbose=true)
@inside sim.flow.σ[I] = WaterLily.curl(3,I,sim.flow.u)*sim.L/sim.U
flood(sim.flow.σ,clims=(-10,10))

# qfield influence
qϕ = -∫2Qϕ!(fpm,sim.flow,0.0,recompute=false)

# we can flood it as the influence is still there
flood(fpm.σ,cmap=:coolwarm)
