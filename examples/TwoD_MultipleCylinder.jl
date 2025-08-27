using WaterLily,StaticArrays,ForcePartition,Plots

function make_circle(;L=32,Re=250,U=1,T=Float32,mem=Array)
    # make a body
    body = AutoBody((x,t)->√sum(abs2, x .- SA[2L,2L] .- 1.5f0) - L/2)
    body += AutoBody((x,t)->√sum(abs2, x .- SA[3.5L,1.5L] .+ 1.5f0) - L/2)

    # generate sim
    Simulation((8L,4L), (U,0), L; ν=U*L/Re, body, T, mem)
end

# run a sim and plot the time evolution
sim = make_circle(L=32); sim_step!(sim, 20)
measure!(sim.flow,sim.body); update!(sim.pois)
# Q-field
Q = copy(sim.flow.σ);
@inside Q[I] = ForcePartition.Qcriterion(I,sim.flow.u).*sim.L/sim.U
# @inside sim.flow.σ[I] = ifelse(abs(sim.flow.σ[I])<0.001,0.0,sim.flow.σ[I])
flood(Q[inside(Q)]; levels=50,lw=0.2)

# force partition
fpm = ForcePartitionMethod(sim)
potential!(fpm,fpm.body;tᵢ=0,axis=1)
active = sim.body.a
passive = sim.body.b
@inside fpm.ϕ[I] =  WaterLily.μ₀(sdf(active,loc(0,I),0),1)*WaterLily.μ₀(sdf(passive,loc(0,I),0),1)*fpm.ϕ[I]
p1 = flood(fpm.ϕ[inside(fpm.ϕ)]; levels=100,clims=(-32,32),lw=0.2,axis=([], false),
           border=:none, title="ϕ₁ | nᵢ⋅∇ϕ=nᵢ on B₁ ∪ B₂")
body_plot!(sim); annotate!([2sim.L,3.5sim.L],[2sim.L,1.5sim.L],["B₁","B₂"])
Qϕ = -2.0.*(fpm.ϕ.*Q); int=round(sum(Qϕ)/sim.L,digits=3)
p12 = flood(Qϕ[inside(Qϕ)],clims=(-4,4),levels=30,lw=0.2,axis=([], false),cfill=:bam,
           border=:non,title="-2∫Qϕ₁dV=$int")
body_plot!(sim); annotate!([2sim.L,3.5sim.L],[2sim.L,1.5sim.L],["B₁","B₂"])
 # new potential only with the active body
potential!(fpm,active;tᵢ=0,axis=1)
@inside fpm.ϕ[I] =  WaterLily.μ₀(sdf(active,loc(0,I),0),1)*WaterLily.μ₀(sdf(passive,loc(0,I),0),1)*fpm.ϕ[I]
ϕ₁_B₁B₂ = copy(fpm.ϕ)
p2 = flood(fpm.ϕ[inside(fpm.ϕ)]; levels=100,clims=(-32,32),lw=0.2,axis=([], false),
           border=:none, title="ϕ₁ | nᵢ⋅∇ϕ=nᵢ on B₁; nᵢ⋅∇ϕ=0 on B₂")
body_plot!(sim); annotate!([2sim.L,3.5sim.L],[2sim.L,1.5sim.L],["B₁","B₂"])
Qϕ = -2.0.*(fpm.ϕ.*Q); int=round(sum(Qϕ)/sim.L,digits=3)
p21 = flood(Qϕ[inside(Qϕ)],clims=(-4,4),levels=30,lw=0.2,axis=([], false),cfill=:bam,
            border=:none,title="-2∫Qϕ₁dV=$int")
body_plot!(sim); annotate!([2sim.L,3.5sim.L],[2sim.L,1.5sim.L],["B₁","B₂"])

# case where the second body is removed entierly
measure!(sim.flow,active); update!(sim.pois)
potential!(fpm,active;tᵢ=0,axis=1)
@inside fpm.ϕ[I] = WaterLily.μ₀(sdf(active,loc(0,I),0),1)*fpm.ϕ[I]
ϕ₁_B₁ = copy(fpm.ϕ)
p3 = flood(fpm.ϕ[inside(fpm.ϕ)]; levels=100,clims=(-32,32),lw=0.2,axis=([], false),
           border=:none, title="ϕ₁ | nᵢ⋅∇ϕ=nᵢ on B₁; B₂ ≡ ∅")
body_plot!(sim); annotate!([2sim.L,3.5sim.L],[2sim.L,1.5sim.L],["B₁","B₂"])
Qϕ = -2.0.*(fpm.ϕ.*Q); int=round(sum(Qϕ)/sim.L,digits=3)
p31 = flood(Qϕ[inside(Qϕ)],clims=(-4,4),levels=30,lw=0.2,axis=([], false),cfill=:bam,
            border=:none,title="-2∫Qϕ₁dV=$int")
body_plot!(sim); annotate!([2sim.L,3.5sim.L],[2sim.L,1.5sim.L],["B₁","B₂"])
# save fig
plot(p1,p12,p2,p21,p3,p31,layout=(3,2),size=(1800,1800))
savefig("assets/potential_influence_various_BCs.png")

# remove the influence inside the bodies
@inside ϕ₁_B₁[I] = WaterLily.μ₀(sdf(passive,loc(0,I),0),1)*ϕ₁_B₁[I]
Qϕ = (ϕ₁_B₁B₂.-ϕ₁_B₁); int=round(sum(Qϕ[inside(Qϕ)])/sim.L,digits=3)
p1=flood(Qϕ[inside(Qϕ)],levels=30,lw=0.2,axis=([],false),
         title="ϕ₁(B₁ ∩ B₂)-ϕ₁(B₁)")
body_plot!(sim); annotate!([2sim.L,3.5sim.L],[2sim.L,1.5sim.L],["B₁","B₂"])
Qϕ = -2Qϕ.*Q; int=round(sum(Qϕ[inside(Qϕ)])/sim.L,digits=3)
p2=flood(Qϕ[inside(Qϕ)],clims=(-4,4),levels=30,lw=0.2,axis=([],false),cfill=:bam,
         title="-2∫Q[ϕ₁(B₁ ∩ B₂)-ϕ₁(B₁)]dV=$int")
body_plot!(sim); annotate!([2sim.L,3.5sim.L],[2sim.L,1.5sim.L],["B₁","B₂"])
plot(p1,p2,layout=(2,1),size=(1280,1280))
savefig("assets/difference_potentials_homogeneousNeumann_nobc.png")