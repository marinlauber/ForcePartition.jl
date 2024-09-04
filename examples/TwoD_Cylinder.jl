using WaterLily
using Plots
using StaticArrays
using ForcePartition
using ParametricBodies
using BiotSavartBCs
include("TwoD_plots.jl")

function make_circle(;L=32,Re=250,U=1)
    radius, center = L/2, 4L
    center = SA[center,center]

    # make a body
    circle = AutoBody((x,t)->√sum(abs2, x .- center .- 1.5) - radius)

    # generate sim
    Simulation((16*L,8*L), (U,0), radius; ν=U*radius/Re, body=circle)
end

# make simulation
sim = make_circle(L=32)

# contruct FPM and compute potential
fpm = ForcePartitionMethod(sim)
potential!(fpm;tᵢ=0,f=:force,axis=1,itmx=32)

# evolve the simulation
sim_gif!(sim;duration=20,step=0.2,verbose=true,
         clims=(-5,5),plotbody=true,remeasure=false)
@inside sim.flow.σ[I] = WaterLily.curl(3,I,sim.flow.u)*sim.L/sim.U
flood(sim.flow.σ,clims=(-5,5))

# plot to check
R = inside(fpm.ϕ)
flood(fpm.ϕ[R],levels=11,filled=true)
body_plot!(sim)
xlims!(100,200); ylims!(100,160)

# gradient should be similar to potential flow solution, except the frame of refernce we are in
R = inside(fpm.ϕ)
@inside sim.flow.σ[I] = 1-WaterLily.∂(1,I,fpm.ϕ)
@inside sim.flow.σ[I] = ifelse(sdf(sim.body,loc(0,I),0)<0.01,NaN,sim.flow.σ[I])
p1 = flood(sim.flow.σ[R],clims=(0,2),levels=11,filled=true,frame=:none,title="U-velocity")
body_plot!(sim); xlims!(100,200); ylims!(100,160)

@inside sim.flow.σ[I] = -WaterLily.∂(2,I,fpm.ϕ) 
@inside sim.flow.σ[I] = ifelse(sdf(sim.body,loc(0,I),0)<0.01,NaN,sim.flow.σ[I])
p2 = flood(sim.flow.σ[R],clims=(-1,1),levels=11,filled=true,frame=:none,title="V-velocity")
body_plot!(sim); xlims!(100,200); ylims!(100,160)

plot(p1,p2, layout = @layout [a b])
# savefig("potential_flow.png")

# test q-criterion
@inside sim.flow.σ[I] = ForcePartition.Qcriterion(I,sim.flow.u)*sim.L/sim.U
flood(sim.flow.σ,clims=(0,1/sim.L),cmap=:blues,lw=0)
xlims!(64,300); ylims!(68,192)

# what is the influemce field
@inside sim.flow.σ[I] = fpm.ϕ[I]*ForcePartition.Qcriterion(I,sim.flow.u)
flood(sim.flow.σ[R],clims=(0,1/sim.L),levels=31,cmap=:blues)
xlims!(64,300); ylims!(68,192)

# qfield influence
qϕ = -∫2Qϕ!(fpm,sim.flow,20.0,recompute=false)
# total force
force = -WaterLily.pressure_force(sim)

# run a sim and plot the time evolution
sim = make_circle(L=32)
fpm = ForcePartitionMethod(sim)
potential!(fpm;tᵢ=0,f=:force,axis=1,itmx=32)
t₀,duration,step = 0.,100,0.2
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
    @inside sim.flow.σ[I] = ifelse(abs(sim.flow.σ[I])<0.001,0.0,sim.flow.σ[I])
    flood(sim.flow.σ; levels=20, shift=(-0.5,-0.5),clims=(-5,5))
    body_plot!(sim); plot!(title="tU/L $tᵢ")

    # print time step
    println("tU/L=",round(tᵢ,digits=4),", Δt=",round(sim.flow.Δt[end],digits=3))
end

# # plot results
# plot(cumsum(@views(sim.flow.Δt[1:end-1]))./sim.L,force/2sim.L,label="Total force",
#      xlabel="tU/L",ylabel="2F/ρU²L",ylims=(0,2),xlims=(0,100))
# plot!(cumsum(@views(sim.flow.Δt[1:end-1]))./sim.L,fp/2sim.L,label="vorticity force")
# savefig("force_partition.png")
