using WaterLily,StaticArrays,ForcePartition,Plots

function make_circle(;L=32,Re=250,U=1,T=Float32,mem=Array)
    # make a body
    body = AutoBody((x,t)->√sum(abs2, x .- SA[2L,2L] .- 1.5f0) - L/2)
    body += AutoBody((x,t)->√sum(abs2, x .- SA[3.5L,1.5L] .+ 1.5f0) - L/2)

    # generate sim
    Simulation((8L,4L), (U,0), L; ν=U*L/Re, body, T, mem)
end

# run a sim and plot the time evolution
# sim = make_circle(L=32); sim_step!(sim, 10)
measure!(sim.flow,sim.body); update!(sim.pois)
# Q-field
Q = copy(sim.flow.σ);
@inside Q[I] = ForcePartition.Qcriterion(I,sim.flow.u)
# @inside sim.flow.σ[I] = ifelse(abs(sim.flow.σ[I])<0.001,0.0,sim.flow.σ[I])
flood(Q[inside(Q)]; levels=50,lw=0.2)

# force partition
fpm = ForcePartitionMethod(sim)
potential!(fpm,fpm.body;tᵢ=0,axis=1)
active = sim.body.a
passive = sim.body.b
@inside fpm.ϕ[I] =  WaterLily.μ₀(sdf(active,loc(0,I),0),1)*WaterLily.μ₀(sdf(passive,loc(0,I),0),1)*fpm.ϕ[I]
p1 = flood(fpm.ϕ[inside(fpm.ϕ)]; levels=100,clims=(-32,32),lw=0.2,axis=([], false),
           legend=false,border=:none, title="ϕ₁ | nᵢ⋅∇ϕ=nᵢ on B₁ ∪ B₂")
body_plot!(sim); annotate!([2sim.L,3.5sim.L],[2sim.L,1.5sim.L],["B₁","B₂"])
Qϕ = -2.0.*(fpm.ϕ.*Q); int=round(sum(Qϕ)\sim.L,digits=3)
p12 = flood(Qϕ[inside(Qϕ)],clims=(-1,1),levels=30,lw=0.2,axis=([], false),cfill=:bam,
           legend=false,border=:non,title="-2∫Qϕ₁dV=$int")
body_plot!(sim); annotate!([2sim.L,3.5sim.L],[2sim.L,1.5sim.L],["B₁","B₂"])
 # new potential only with the active body
potential!(fpm,active;tᵢ=0,axis=1)
@inside fpm.ϕ[I] =  WaterLily.μ₀(sdf(active,loc(0,I),0),1)*WaterLily.μ₀(sdf(passive,loc(0,I),0),1)*fpm.ϕ[I]
ϕ_1_neumann = copy(fpm.ϕ)
p2 = flood(fpm.ϕ[inside(fpm.ϕ)]; levels=100,clims=(-32,32),lw=0.2,axis=([], false),
           legend=false,border=:none, title="ϕ₁ | nᵢ⋅∇ϕ=nᵢ on B₁; nᵢ⋅∇ϕ=0 on B₂")
body_plot!(sim); annotate!([2sim.L,3.5sim.L],[2sim.L,1.5sim.L],["B₁","B₂"])
Qϕ = -2.0.*(fpm.ϕ.*Q); int=round(sum(Qϕ)\sim.L,digits=3)
p21 = flood(Qϕ[inside(Qϕ)],clims=(-1,1),levels=30,lw=0.2,axis=([], false),cfill=:bam,
           legend=false,border=:none,title="-2∫Qϕ₁dV=$int")
body_plot!(sim); annotate!([2sim.L,3.5sim.L],[2sim.L,1.5sim.L],["B₁","B₂"])

# case where the second body is removed entierly
measure!(sim.flow,active); update!(sim.pois)
potential!(fpm,active;tᵢ=0,axis=1)
@inside fpm.ϕ[I] =  WaterLily.μ₀(sdf(active,loc(0,I),0),1)*fpm.ϕ[I]
ϕ_1_nob2 = copy(fpm.ϕ)
p3 = flood(fpm.ϕ[inside(fpm.ϕ)]; levels=100,clims=(-32,32),lw=0.2,axis=([], false),
           legend=false,border=:none, title="ϕ₁ | nᵢ⋅∇ϕ=nᵢ on B₁; B₂ ≡ ∅")
Qϕ = -2.0.*(fpm.ϕ.*Q); int=round(sum(Qϕ)\sim.L,digits=3)
p31 = flood(Qϕ[inside(Qϕ)],clims=(-1,1),levels=30,lw=0.2,axis=([], false),cfill=:bam,
            legend=false,border=:none,title="-2∫Qϕ₁dV=$int")
body_plot!(sim); annotate!([2sim.L,3.5sim.L],[2sim.L,1.5sim.L],["B₁","B₂"])
# save fig
plot(p1,p12,p2,p21,p3,p31,layout=(3,2),size=(1800,1800))
savefig("potential_influence_various_BCs.png")

Qϕ = -2Q.*((ϕ_1_neumann.-ϕ_1_nob2)); int=round(sum(Qϕ[inside(Qϕ)])\sim.L,digits=3)
flood(Qϕ[inside(Qϕ)],levels=30,lw=0.2,axis=([],false),cfill=:bam,
      legend=false,border=:none,title="-2∫Qϕ₁dV=$int")
body_plot!(sim); annotate!([2sim.L,3.5sim.L],[2sim.L,1.5sim.L],["B₁","B₂"])
savefig("difference_potentials_homogeneousNeumann_nobc.png")
