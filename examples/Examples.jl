using WaterLily,StaticArrays,ForcePartition

function make_circle(;L=32,Re=250,U=1)
    # parameters
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
potential!(fpm,fpm.body;tᵢ=0,axis=1)

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
@inside sim.flow.σ[I] = 1-ForcePartition.∇ϕ(I,fpm.ϕ)[1]
# @inside sim.flow.σ[I] = 1-WaterLily.∂(1,I,fpm.ϕ)
@inside sim.flow.σ[I] = ifelse(sdf(sim.body,loc(0,I),0)<0.01,NaN,sim.flow.σ[I])
p1 = flood(sim.flow.σ[R],clims=(0,2),levels=11,filled=true,frame=:none,title="U-velocity")
body_plot!(sim); xlims!(100,200); ylims!(100,160)

@inside sim.flow.σ[I] = -ForcePartition.∇ϕ(I,fpm.ϕ)[2]
# @inside sim.flow.σ[I] = -WaterLily.∂(2,I,fpm.ϕ)
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
qϕ = -∫2QϕdV!(fpm,sim.flow,20.0,recompute=false)
# total force
force = -WaterLily.pressure_force(sim)

function make_sphere(;L=32,Re=250,U=1)
    radius, center = L/2, SA[2L,2L,2L]
    # make a body
    sphere = AutoBody((x,t)->√sum(abs2, x .- center .- 1.5) - radius)

    # generate sim
    Simulation((8L,4L,4L), (U,0,0), radius; ν=U*radius/Re, body=sphere)
end

# make simulation
sim = make_sphere(L=32)

# contruct FPM and compute potential
fpm = ForcePartitionMethod(sim)
potential!(fpm,fpm.body;tᵢ=0,f=:force,axis=1,itmx=32)

# evolve the simulation
sim_step!(sim,1)
@inside sim.flow.σ[I] = WaterLily.curl(2,I,sim.flow.u)*sim.L/sim.U
flood(sim.flow.σ[:,:,64],clims=(-5,5))

# check the influence field
contourf(fpm.ϕ[:,:,64]',aspect_ratio=:equal,cmap=:seismic)
contourf(fpm.ϕ[:,64,:]',aspect_ratio=:equal,cmap=:seismic)

# gradient should be similar to potential flow solution, except the frame of refernce we are in
@inside sim.flow.σ[I] = 1-WaterLily.∂(1,I,fpm.ϕ)
@inside sim.flow.σ[I] = ifelse(sdf(sim.body,loc(0,I),0)<0.01,NaN,sim.flow.σ[I])
p1 = contourf(sim.flow.σ[:,:,64]',clims=(0,2),aspect_ratio=:equal,levels=11,filled=true,frame=:none,title="U-velocity")

@inside sim.flow.σ[I] = -WaterLily.∂(2,I,fpm.ϕ)
@inside sim.flow.σ[I] = ifelse(sdf(sim.body,loc(0,I),0)<0.01,NaN,sim.flow.σ[I])
p2 = contourf(sim.flow.σ[:,:,64]',clims=(-1,1),aspect_ratio=:equal,levels=11,filled=true,frame=:none,title="V-velocity")
plot(p1,p2, layout = @layout [a b])

# qfield influence
qϕ = -∫2Qϕ!(fpm,sim.flow,recompute=false)
# total force
force = -WaterLily.pressure_force(sim)