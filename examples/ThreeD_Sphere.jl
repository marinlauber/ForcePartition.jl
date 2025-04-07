using WaterLily,StaticArrays,ForcePartition,Plots

function make_sphere(;L=32,Re=250,U=1)
    radius, center = L/2, SA[2L,2L,2L]
    # make a body
    sphere = AutoBody((x,t)->√sum(abs2, x .- center .- 1.5) - radius)

    # generate sim
    Simulation((8L,4L,4L), (U,0,0), radius; ν=U*radius/Re, body=sphere)
end

# run a sim and plot the time evolution
sim = make_sphere(L=32)
fpm = ForcePartitionMethod(sim)
potential!(fpm,fpm.body;tᵢ=0,axis=1,itmx=32)
t₀,duration,step = 0.,10,0.2
force, fp = [],[]

@time for tᵢ in range(t₀,t₀+duration;step)
    # update until time tᵢ in the background
    while sim_time(sim) < tᵢ
        # update flow
        mom_step!(sim.flow,sim.pois)
        # pressure force
        push!(force,-2WaterLily.pressure_force(sim)[1])
        push!(fp,-∫2QϕdV!(fpm,sim.flow,recompute=false))
    end

    # print time step
    println("tU/L=",round(tᵢ,digits=4),", Δt=",round(sim.flow.Δt[end],digits=3))
end

# plot results
plot(cumsum(@views(sim.flow.Δt[1:end-1]))./sim.L,force/sim.L^2,label="Total force",
     xlabel="tU/L",ylabel="2F/ρU²L",ylims=(0,2),xlims=(0,10))
plot!(cumsum(@views(sim.flow.Δt[1:end-1]))./sim.L,fp/sim.L^2,label="vorticity force")
savefig("force_partition_sphere.png")