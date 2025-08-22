using WaterLily,StaticArrays,Plots

function make_circle(;L=32,Re=250,U=1,T=Float32,mem=Array)
    # parameters
    radius, center = T(L/2), T(4L)
    center = SA{T}[center,center]

    # make a body
    circle = AutoBody((x,t)->√sum(abs2, x .- center .- 1.5f0) - radius)

    # generate sim
    Simulation((16*L,8*L), (U,0), radius; ν=U*radius/Re, body=circle, T, mem)
end

# run a sim and plot the time evolution
sim = make_circle(L=32)
t₀,duration,step = 0.,100,0.2

@time @gif for tᵢ in range(t₀,t₀+duration;step)
    # update until time tᵢ in the background
    while sim_time(sim) < tᵢ
        # update flow
        mom_step!(sim.flow,sim.pois)
    end

    # plot vorticity
    @inside sim.flow.σ[I] = WaterLily.curl(3,I,sim.flow.u)*sim.L/sim.U
    @inside sim.flow.σ[I] = ifelse(abs(sim.flow.σ[I])<0.001,0.0,sim.flow.σ[I])
    flood(sim.flow.σ; levels=20, shift=(-0.5,-0.5),clims=(-5,5))
    body_plot!(sim); plot!(title="tU/L $tᵢ")

    # print time step
    println("tU/L=",round(tᵢ,digits=4),", Δt=",round(sim.flow.Δt[end],digits=3))
end