using WaterLily,StaticArrays,ForcePartition,ParametricBodies,Plots

function make_airfoil(;L=32,Re=1000,St=0.25,αₘ=25,U=1,n=16,m=8,T=Float32,mem=Array)
    # Map from simulation coordinate x to surface coordinate ξ
    center,pivot = SA[4L,0.5f0m*L],SA[0.5f0L,0]
    θ₀ = T(αₘ*π/180); ω=T(2π*St*U/L)
    function map(x,t)
        θ = θ₀*sin(ω*t); R = SA[cos(θ) -sin(θ); sin(θ) cos(θ)]
        ξ = R*(x-center-pivot)+pivot # move to origin and align with x-axis
        return SA[ξ[1],abs(ξ[2])]    # reflect to positive y
    end

    # Define foil using NACA0012 profile equation: https://tinyurl.com/NACA00xx
    NACA(s) = 0.6f0*(0.2969f0s-0.126f0s^2-0.3516f0s^4+0.2843f0s^6-0.1036f0s^8)
    foil(s,t) = L*SA[(1-s)^2,NACA(1-s)]
    body = HashedBody(foil,(0,1);map,T,mem)

    # make the sim
    Simulation((n*L,m*L),(U,0),L;ν=U*L/Re,body,T,mem)
end

# run a sim and plot the time evolution
sim = make_airfoil(L=32)
fpm = ForcePartitionMethod(sim;axis=6) # z-moment
potential!(fpm,fpm.body;tᵢ=0)
t₀,duration,step = 0.,2.0,0.2
force, fp = [],[]

@time @gif for tᵢ in range(t₀,t₀+duration;step)
    # update until time tᵢ in the background
    while sim_time(sim) < tᵢ
        # update flow
        measure!(sim); mom_step!(sim.flow, sim.pois)
        # pressure force
        push!(force,-2WaterLily.pressure_force(sim)[1])
        push!(fp,-∫2QϕdV!(fpm,sim.flow,recompute=true))
    end

    # plot vorticity
    # @inside sim.flow.σ[I] = WaterLily.curl(3,I,sim.flow.u)*sim.L/sim.U
    @inside sim.flow.σ[I] = fpm.ϕ[I]
    flood(sim.flow.σ; shift=(-0.5,-0.5),clims=(-5,5))
    body_plot!(sim); plot!(title="tU/L $tᵢ")

    # print time step
    println("tU/L=",round(tᵢ,digits=4),", Δt=",round(sim.flow.Δt[end],digits=3))
end

# plot results
plot(cumsum(@views(sim.flow.Δt[1:end-1]))./sim.L,force/2sim.L,label="Total force",
     xlabel="tU/L",ylabel="2F/ρU²L",ylims=(0,2),xlims=(0,100))
plot!(cumsum(@views(sim.flow.Δt[1:end-1]))./sim.L,fp/2sim.L,label="vorticity force")
savefig("force_partition_airfoil.png")