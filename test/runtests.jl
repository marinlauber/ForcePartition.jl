using ForcePartition
using WaterLily
using Test
using StaticArrays

arrays = [Array]
# CUDA.functional() && push!(arrays, CuArray)

# potential around a cylinder
ϕ_cylinder(center,R) = (y-x.-center; r=√sum(abs2,y); θ=atan(y[2]/y[1]); ϕ = -R^3*cos(θ)/2r^2; return ϕ)

# circle sim
function circle(;L=32,m=10,Re=250,U=1,T=Float32,mem=Array)
    # body
    body = AutoBody((x,t)->√sum(abs2, x.-m/2*L)-L÷2)
    # generate sim
    Simulation((m*L,m*L), (U,0), L; ν=U*radius/Re, body, T, mem)
end

@testset "ForcePartitionMethod.jl" begin
    #
    sim = circle(L=32,T=Float32,mem=Array)
    fpm = ForcePartitionMethod(sim)
    potential!(fpm,fpm.body;tᵢ=0,axis=1)
    R = inside(fpm.ϕ)
    @test maximum(abs.(fpm.ϕ[R] .- ϕ_cylinder(SA[sim.L*sim.m/2,sim.L*sim.m/2],sim.L/2)[R])) < 1e-2
end

function parabolic_u(i,x,t,L)
    i ≠ 1 && return 0
    r = √sum(abs2,SA[x[2],x[3]].-L/2)
    return r<L/2-3/2 ? 2*(1-r^2/(L/2-3/2)^2) : 0 # remove radius and add stenosis (and the ghost)
end

@testset "AverageFlow.jl" begin
    # sizes
    N = 64
    for f ∈ arrays
        u = zeros(N+2,N+2,N+2,3) |> f;
        apply!((i,x)->parabolic_u(i,x,0.,N),u)
        x,v = ForcePartition.azimuthal_avrg(@views(u[10,:,:,1]))
        # error
        para(r;L=32) = r<L/2-3/2 ? 2*(1-r^2/(L/2-3/2)^2) : 0
        ϵ(r,u;N=32) = √sum(abs2,[para(rᵢ;L=N)-uᵢ for (rᵢ,uᵢ) in zip(r,u)])/length(u)
        @test ϵ(x,v;N) < 1e-2
    end
end
