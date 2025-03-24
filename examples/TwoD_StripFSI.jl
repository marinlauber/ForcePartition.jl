using WaterLily,StaticArrays,ForcePartition,Splines,Plots

include("/home/marin/Workspace/Tutorials-WaterLily/src/TwoD_plots.jl")

function circle(D,dims,Ubc;Re=200,U=1,mem=Array)
    body = AutoBody((x,t)->√sum(abs2,SA[x[1],x[2]].-dims[2]÷2)-D÷2)
    Simulation(dims, Ubc, D; body, ν=U*D/Re, mem)
end

function step!(sim::NTuple{N,Simulation};Δt_min=Inf) where N
    for i in 1:N
        mom_step!(sim[i].flow,sim[i].pois)
        Δt_min = min(Δt_min,sim[i].flow.Δt[end])
    end
    # keep time in sync
    foreach(i->sim[i].flow.Δt[end]=Δt_min,N)
end

# make flow a collection of strips
Ns = 3
R = 16
sim = ntuple(i->circle(2R,(20R,8R),(1,0);Re=10_000,U=1,mem=Array), Ns)

# initialize and run
@gif for t in range(0,20.0;step=0.05)#1:6
    while sim_time(sim[1])<t #sim_step!(sim,t)
        step!(sim) # strip update
    end
    @inside sim[1].flow.σ[I] = WaterLily.curl(3,I,sim[1].flow.u)*sim[1].L/sim[1].U
    p1 = flood(sim[1].flow.σ,shift=(-2,-1.5),clims=(-5,5))
    @inside sim[2].flow.σ[I] = WaterLily.curl(3,I,sim[2].flow.u)*sim[2].L/sim[2].U
    p2 = flood(sim[2].flow.σ,shift=(-2,-1.5),clims=(-5,5))
    @inside sim[3].flow.σ[I] = WaterLily.curl(3,I,sim[3].flow.u)*sim[3].L/sim[3].U
    p3 = flood(sim[3].flow.σ,shift=(-2,-1.5),clims=(-5,5))
    plot(p1,p2,p3,layout=@layout [a b c])
    @show t
end