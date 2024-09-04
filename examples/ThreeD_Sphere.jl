using WaterLily
using Plots
using StaticArrays
using ForcePartition
using ParametricBodies
include("TwoD_plots.jl")

function make_sphere(;L=32,Re=250,U=1)
    radius, center = L/2, SA[2L,2L,2L]
    # make a body
    sphere = AutoBody((x,t)->√sum(abs2, x .- center .- 1.5) - radius)

    # generate sim
    Simulation((8L,4L,4L), (U,0,0), radius; ν=U*radius/Re, body=sphere)
end

# # make simulation
# sim = make_sphere(L=32)

# # contruct FPM and compute potential
# fpm = ForcePartitionMethod(sim)
# potential!(fpm;tᵢ=0,f=:force,axis=1,itmx=32)

# # evolve the simulation
# sim_step!(sim,1)
# @inside sim.flow.σ[I] = WaterLily.curl(2,I,sim.flow.u)*sim.L/sim.U
# flood(sim.flow.σ[:,:,64],clims=(-5,5))

# # check the influence field
# contourf(fpm.ϕ[:,:,64]',aspect_ratio=:equal,cmap=:seismic)
# contourf(fpm.ϕ[:,64,:]',aspect_ratio=:equal,cmap=:seismic)

# # gradient should be similar to potential flow solution, except the frame of refernce we are in
# @inside sim.flow.σ[I] = 1-WaterLily.∂(1,I,fpm.ϕ)
# @inside sim.flow.σ[I] = ifelse(sdf(sim.body,loc(0,I),0)<0.01,NaN,sim.flow.σ[I])
# p1 = contourf(sim.flow.σ[:,:,64]',clims=(0,2),aspect_ratio=:equal,levels=11,filled=true,frame=:none,title="U-velocity")

# @inside sim.flow.σ[I] = -WaterLily.∂(2,I,fpm.ϕ) 
# @inside sim.flow.σ[I] = ifelse(sdf(sim.body,loc(0,I),0)<0.01,NaN,sim.flow.σ[I])
# p2 = contourf(sim.flow.σ[:,:,64]',clims=(-1,1),aspect_ratio=:equal,levels=11,filled=true,frame=:none,title="V-velocity")
# plot(p1,p2, layout = @layout [a b])

# # qfield influence
# qϕ = -∫2Qϕ!(fpm,sim.flow,recompute=false)
# # total force
# force = -WaterLily.pressure_force(sim)

# run a sim and plot the time evolution
sim = make_sphere(L=32)
fpm = ForcePartitionMethod(sim)
potential!(fpm;tᵢ=0,f=:force,axis=1,itmx=32)
t₀,duration,step = 0.,10,0.2
force, fp = [],[]

@time for tᵢ in range(t₀,t₀+duration;step)
    # update until time tᵢ in the background
    while sim_time(sim) < tᵢ
        # update flow
        mom_step!(sim.flow,sim.pois)
        # pressure force
        push!(force,-2WaterLily.pressure_force(sim)[1])
        push!(fp,-∫2Qϕ!(fpm,sim.flow,recompute=false))
    end

    # print time step
    println("tU/L=",round(tᵢ,digits=4),", Δt=",round(sim.flow.Δt[end],digits=3))
end

# # plot results
# plot(cumsum(@views(sim.flow.Δt[1:end-1]))./sim.L,force/sim.L^2,label="Total force",
#      xlabel="tU/L",ylabel="2F/ρU²L",ylims=(0,2),xlims=(0,10))
# plot!(cumsum(@views(sim.flow.Δt[1:end-1]))./sim.L,fp/sim.L^2,label="vorticity force")
# savefig("force_partition_sphere.png")