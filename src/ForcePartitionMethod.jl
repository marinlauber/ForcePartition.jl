using WaterLily
using StaticArrays
using LinearAlgebra: tr

Base.copy(b::AutoBody) = AutoBody(b.sdf,b.map)
struct ForcePartitionMethod{T,Sf<:AbstractArray{T},Vf<:AbstractArray{T}}
    ϕ :: Sf
    σ :: Sf
    f :: Vf
    pois::AbstractPoisson
    body::AbstractBody
    function ForcePartitionMethod(sim::Simulation)
        T = eltype(sim.flow.u)
        ϕ = copy(sim.flow.p)
        σ = copy(sim.flow.σ)
        f = copy(sim.flow.f)
        body = copy(sim.body)
        new{T,typeof(ϕ),typeof(f)}(ϕ,σ,f,Poisson(ϕ,sim.flow.μ₀,σ),body)
    end
end

"""
    Compute the the influence potential of the body at time tᵢ 
"""
# #@TODO check that updating the body updates the poisson coefficients of the FPM
function potential!(FPM::ForcePartitionMethod,tᵢ=0;axis=1)
    # the boundary value is the surface normal
    function normal(x,i)
        dᵢ,nᵢ,_ = measure(FPM.body,x,tᵢ)
        μ₀ = WaterLily.μ₀(dᵢ,1)
        (1.0-μ₀)*nᵢ[i] # i-component of the normal
    end
    # generate source term
    apply!(x->normal(x,axis),FPM.pois.z); BC!(FPM.pois.z)
    # solver for potential
    solver!(FPM.pois); BC!(FPM.ϕ)
end

function Qcriterion(I::CartesianIndex{2},u)
    J = @SMatrix [WaterLily.∂(i,j,I,u) for i ∈ 1:2, j ∈ 1:2]
    S,Ω = (J+J')/2,(J-J')/2
    ## -0.5*sum(eigvals(S^2+Ω^2)) # this is also possible, but 2x slower
    0.5*(√(tr(Ω*Ω'))^2-√(tr(S*S'))^2)
end

function ∫2Qϕ!(FPM::ForcePartitionMethod,a::Flow,tᵢ=0;axis=1,recompute=true)
    # get potential
    recompute && potential!(FPM,tᵢ;axis=axis)
    # compute the influence of the Q field
    @WaterLily.loop FPM.σ[I] = FPM.ϕ[I]*Qcriterion(I,a.u) over I ∈ inside(FPM.σ)
    # return the integral over the doman
    sum(@inbounds(FPM.σ[inside(FPM.σ)]))
end