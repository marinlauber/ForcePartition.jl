using ForcePartition
using WaterLily
using Test
using StaticArrays

arrays = [Array]
# CUDA.functional() && push!(arrays, CuArray)

@testset "ForcePartition.jl" begin
    # Write your tests here.
    @test true
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
