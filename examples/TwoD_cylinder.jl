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
    body = ParametricBody(foil,(0,1);map)

    # make the sim
    Simulation((n*L,m*L),(U,0),L;ν=U*L/Re,body,T,mem)
end


# make simulation
# sim = make_circle(L=32)
sim = make_airfoil(L=32)

# contruct fpm and compute potential
fpm = ForcePartitionMethod(sim)
potential!(fpm;tᵢ=0,f=:moment,axis=1,itmx=1e6)

# plot to check
R = inside(fpm.ϕ)
flood(fpm.ϕ[R],levels=11,filled=true)
body_plot!(sim)
xlims!(100,200); ylims!(100,150)

# evolve the simulation
sim_gif!(sim;duration=20,step=0.2,verbose=true,
         clims=(-5,5),plotbody=true,remeasure=true)
@inside sim.flow.σ[I] = WaterLily.curl(3,I,sim.flow.u)*sim.L/sim.U
flood(sim.flow.σ,clims=(-10,10))

# plot to check
R = inside(fpm.ϕ)
flood(fpm.ϕ[R],clims=(-sim.L,sim.L),levels=11,filled=true)
body_plot!(sim)


# qfield influence
potential!(fpm;x₀=SA[4.5sim.L,4sim.L],tᵢ=20,f=:force,axis=1)
qϕ = -∫2Qϕ!(fpm,sim.flow,20.0,recompute=false)

force = -WaterLily.∮nds(sim.flow.p,sim.flow.f,sim.body,WaterLily.time(sim))

# we can flood it as the influence is still there
flood(fpm.ϕ,clims=(-100,0),cmap=:RdBu)
flood(fpm.σ,cmap=:RdBu,clims=(-0.25,0.25))
# plot!([4.5sim.L],[4sim.L],marker=:circle,markersize=5,markercolor=:black)
body_plot!(sim)

# plot Q-criterion on its own
@inside sim.flow.σ[I] = ForcePartition.Qcriterion(I,sim.flow.u)*sim.L/sim.U
flood(clamp.(sim.flow.σ,-Inf,Inf),cmap=:Reds)
body_plot!(sim)

#plot the source term
fpm = ForcePartitionMethod(sim,MultiLevelPoisson)
potential!(fpm;tᵢ=1,f=:force,axis=1)

# plot to check
R = inside(fpm.ϕ)
flood(fpm.ϕ[R],clims=(-sim.L,sim.L),levels=11,filled=true)
body_plot!(sim)

# apply!(x->ForcePartition.force(fpm.body,x,SA[1.5sim.L,2sim.L],1,0.1),fpm.pois.z)
# apply!(x->ForcePartition.moment(fpm.body,x,SA[1.5sim.L,2sim.L],1,5),fpm.pois.z)
# BC!(fpm.pois.z)
# flood(fpm.pois.z,levels=10,shift=(-1.,-1.))
# WaterLily.measure_sdf!(sim.flow.σ,sim.body,0.1)
# contour!(sim.flow.σ[inside(sim.flow.σ)]'|>Array;levels=[0],lines=:black)
# xlims!(30,64);ylims!(60,70)
# solver!(fpm.pois,itmx=32)
# flood(fpm.pois.x,clims=(-1,1))
