using WaterLily,StaticArrays,ForcePartition,ParametricBodies,Plots

function make_airfoil(;L=32,Re=1000,St=0.5,αₘ=25,U=1,T=Float32,mem=Array)
    # Map from simulation coordinate x to surface coordinate ξ
    center,pivot = SA[3L,3L],SA[0.25f0L,0]
    θ₀ = T(αₘ*π/180); ω=T(2π*St*U/L)
    function map(x,t)
        θ = θ₀*sin(ω*t); R = SA[cos(θ) -sin(θ); sin(θ) cos(θ)]
        ξ = R*(x-center-pivot)+pivot # move to origin and align with x-axis
        return SA[ξ[1],abs(ξ[2])]    # reflect to positive y
    end

    # Define foil using NACA0012 profile equation: https://tinyurl.com/NACA00xx
    NACA(s) = 1.2f0*(0.2969f0s-0.126f0s^2-0.3516f0s^4+0.2843f0s^6-0.1036f0s^8)
    foil(s,t) = L*SA[(1-s)^2,NACA(1-s)]
    body = HashedBody(foil,(0,1);map,T,mem)

    # make the sim an return center of rotation
    Simulation((10L,6L),(U,0),L;ν=U*L/Re,body,T,mem),center+pivot
end

# run a sim and plot the time evolution
# using CUDA
sim,xC = make_airfoil(L=64,mem=Array)
fpm = ForcePartitionMethod(sim;axis=6) # z-moment
potential!(fpm;x₀=xC)

# sim setup
t₀,duration,step = 0.,20.0,0.05
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
        push!(fp,-∫2QϕdV!(fpm,sim.flow,x₀=xC,recompute=true))
        push!(fv, 0.0)
    end
    # plot -2Qψ
    potential!(fpm;x₀=xC,tᵢ=WaterLily.time(sim))
    @inside Qϕ[I] = -2.0fpm.ϕ[I]*ForcePartition.Qcriterion(I,sim.flow.u)
    flood(Qϕ[inside(Qϕ)];clims=(-1,1),levels=20,axis=([],false),cfill=:bam,border=:none)
    body_plot!(sim); scatter!([xC[1]],[xC[2]],label="xₔ",ms=3,mc=:black,title="-2Qψ₃")

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