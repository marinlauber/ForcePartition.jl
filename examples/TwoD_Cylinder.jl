using WaterLily,StaticArrays,ForcePartition,Plots

function make_circle(;L=32,Re=550,U=1,T=Float32,mem=Array)
    # make a body
    body = AutoBody((x,t)->√sum(abs2,x)-L÷2.f0,
                    (x,t)->x-SA[3.f0L,min(2.8f0L+t/4.f0,3.f0L)]) # reflect to positive y

    # generate sim
    Simulation((10L,6L),(U,0),L;ν=U*L/Re,body,T,mem)
end

# run a sim and plot the time evolution
using CUDA
sim = make_circle(L=64,mem=Array)
fpm = ForcePartitionMethod(sim)
potential!(fpm,fpm.body;tᵢ=0,axis=1)
t₀,duration,step = 0.,50,0.05
force,forcev,fp,fv = [],[],[],[]
Qϕ = sim.flow.σ # pointer

@time @gif for tᵢ in range(t₀,t₀+duration;step)
    # update until time tᵢ in the background
    while sim_time(sim) < tᵢ
        # update flow
        sim_step!(sim;remeasure=true)
        # pressure force
        push!(force,-2WaterLily.pressure_force(sim)[1])
        push!(forcev,-2WaterLily.viscous_force(sim)[1])
        push!(fp,-∫2QϕdV!(fpm,sim.flow,recompute=true))
        push!(fv, ∮ReωdS!(fpm,sim.flow,recompute=false))
    end

    # plot -2Qϕ
    @inside Qϕ[I] = -2.0fpm.ϕ[I]*ForcePartition.Qcriterion(I,sim.flow.u)
    flood(Qϕ[inside(Qϕ)];clims=(-.1,.1),levels=20,axis=([],false),cfill=:bam,border=:none)
    body_plot!(sim)

    # print time step
    println("tU/L=",round(tᵢ,digits=4),", Δt=",round(sim.flow.Δt[end],digits=3))
end

# plot results
time = cumsum(@views(sim.flow.Δt[1:end-1]))./sim.L
p1=plot(time,[forcev/2sim.L,-fv/2sim.L],
        label=["viscous force" "viscous force (partition)"],
        xlabel="tU/L",ylabel="2F/ρU²L",ylims=(0,0.1),xlims=(20,50),dpi=600)
p2=plot(time,[force/2sim.L,fp/2sim.L,(fp.-fv)/2sim.L],
        label=["pressure force" "vorticity force (partition)" "total partition force"],
        xlabel="tU/L",ylabel="2F/ρU²L",ylims=(0,1),xlims=(20,50),dpi=600)
plot(p1,p2,layout=(1,2),size=(800,400))
savefig("assets/force_partition.png")
