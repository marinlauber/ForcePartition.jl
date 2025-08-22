using WaterLily,StaticArrays,ForcePartition,Plots

function make_circle(;L=32,Re=250,U=1,T=Float32,mem=Array)
    # make a body
    body = AutoBody((x,t)->√sum(abs2, x .- SA[2L,2L] .- 1.5f0) - L/2)
    body += AutoBody((x,t)->√sum(abs2, x .- SA[3.5L,1.5L] .+ 1.5f0) - L/2)

    # generate sim
    Simulation((8L,4L), (U,0), L; ν=U*L/Re, body, T, mem)
end

# run a sim and plot the time evolution
sim = make_circle(L=32)
sim_step!(sim, 10)
# flood
@inside sim.flow.σ[I] = WaterLily.curl(3,I,sim.flow.u)*sim.L/sim.U
@inside sim.flow.σ[I] = ifelse(abs(sim.flow.σ[I])<0.001,0.0,sim.flow.σ[I])
flood(sim.flow.σ; levels=20, shift=(-0.5,-0.5),clims=(-5,5))
# force partition
fpm = ForcePartitionMethod(sim)
potential!(fpm,fpm.body;tᵢ=0,axis=1)
flood(fpm.ϕ; levels=20)












# t₀,duration,step = 0.,10,0.2
# force,forcev,fp,fa,fv = [],[],[],[],[]

# @time @gif for tᵢ in range(t₀,t₀+duration;step)
#     # update until time tᵢ in the background
#     while sim_time(sim) < tᵢ
#         # update flow
#         mom_step!(sim.flow,sim.pois)
#         # pressure force
#         push!(force,-2WaterLily.pressure_force(sim)[1])
#         push!(forcev,-2WaterLily.viscous_force(sim)[1])
#         push!(fa,-∮UϕdS!(fpm,sim.flow,recompute=false))
#         push!(fp,-∫2QϕdV!(fpm,sim.flow,recompute=false))
#         push!(fv, ∮ReωdS!(fpm,sim.flow,recompute=false))
#     end

#     # plot vorticity
#     @inside sim.flow.σ[I] = WaterLily.curl(3,I,sim.flow.u)*sim.L/sim.U
#     @inside sim.flow.σ[I] = ifelse(abs(sim.flow.σ[I])<0.001,0.0,sim.flow.σ[I])
#     flood(sim.flow.σ; levels=20, shift=(-0.5,-0.5),clims=(-5,5))
#     body_plot!(sim); plot!(title="tU/L $tᵢ")

#     # print time step
#     println("tU/L=",round(tᵢ,digits=4),", Δt=",round(sim.flow.Δt[end],digits=3))
# end

# # plot results
# time = cumsum(@views(sim.flow.Δt[1:end-1]))./sim.L
# p1=plot(time,[forcev/2sim.L,-fv/2sim.L],
#         label=["viscous force" "viscous force (partition)"],
#         ylabel="2F/ρU²L",ylims=(0,0.1),xlims=(0,100),dpi=600)
# p2=plot(time,[force/2sim.L,fp/2sim.L,fa/2sim.L,(fp.+fa.-fv)/2sim.L],
#         label=["pressure force" "vorticity force (partition)" "added-mass force (partition)" "total partition force"],
#         xlabel="tU/L",ylabel="2F/ρU²L",ylims=(0,2),xlims=(0,100),dpi=600)
# plot(p1,p2,layout=(2,1),size=(800,600))
# savefig("force_partition.png")
