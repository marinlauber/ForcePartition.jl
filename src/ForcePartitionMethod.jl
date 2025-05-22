using WaterLily
using StaticArrays
using LinearAlgebra: tr,dot

"""
    ForcePartitionMethod(sim::WaterLily.Simulation)
    
    A struct to compute the influence of the flow field on a body.
    - ϕ: potential field (allocates a new field)
    - σ: scalar field (points to `sim.flow.σ` to avoid allocating)
    - f: vector field (points to `sim.flow.f` to avoid allocating)
    - pois: poisson solver (points to `sim.pois` to avoid allocating)
    - body: body (points to `sim.body` to avoid allocating)
"""
struct ForcePartitionMethod{A,T,Sf<:AbstractArray{T},Vf<:AbstractArray{T}}
    ϕ :: Sf
    σ :: Sf
    f :: Vf
    pois :: AbstractPoisson
    body :: AbstractBody
    function ForcePartitionMethod(sim::AbstractSimulation;axis=1)
        ϕ,σ,f = copy(sim.flow.p),sim.flow.σ,sim.flow.f # σ and f point to the flow fields
        new{axis,eltype(sim.flow.u),typeof(ϕ),typeof(f)}(ϕ,σ,f,sim.pois,sim.body) # we point to the poisson solver of the flow
    end
end

"""
    potential!(FPM::ForcePartitionMethod;axis,tᵢ)

    Compute the the influence potential of the body at time tᵢ for either 
    the `:force` or the `:moment`.
"""
function potential!(FPM::ForcePartitionMethod{A,T},body;x₀=0,axis=nothing,tᵢ=0) where {A,T}
    @assert !(ndims(FPM.ϕ)==2 && axis ∈ [3,4,5]) "Only 1, 2 and 6 are valid axis for 2D sims"
    # first we need to save the pressure field
    @inside FPM.σ[I] = FPM.pois.x[I]
    # generate source term
    isnothing(axis) && (axis=A) # if we provide an axis, we use it
    @inside FPM.pois.z[I] = source(body,loc(0,I,T),x₀,axis,tᵢ)
    # solver for potential
    solver!(FPM.pois); pop!(FPM.pois.n) # keep the tol the same as the pressure and don't write the iterations
    # copy to the FPM
    @inside FPM.ϕ[I] = FPM.pois.x[I]
    @inside FPM.pois.x[I] = FPM.σ[I] # put back pressure field
end
potential!(FPM;kwargs...) = potential!(FPM,FPM.body;kwargs...) # for backward compatibility

"""
    source term for the potential, either force or moment depending on the axis
"""
source(body,x,x₀,axis,tᵢ) = axis ≤ 3 ? force(body,x,axis,tᵢ) : moment(body,x,x₀,(axis-1)%3+1,tᵢ)
"""
    force(body,x₀,i,tᵢ)
    
    Compute the the source term for the force influence potential
    of the body at time tᵢ.
"""
function force(body,x,i,tᵢ)
    dᵢ,nᵢ,_ = measure(body,x,tᵢ)
    WaterLily.kern(clamp(dᵢ,-1,1))*nᵢ[i] # i-component of the normal
end
"""
    moment(body,x,x₀,i,tᵢ)

    Compute the the source term for the moment influence potential
    of the body at time tᵢ aroudn point x₀.
"""
function moment(body,x,x₀,i,tᵢ)
    dᵢ,nᵢ,_ = measure(body,x,tᵢ)
    WaterLily.kern(clamp(dᵢ,-1,1))*cross(i,(x.-x₀),nᵢ) # the normal moment
end
cross(i,a::SVector{2},b::SVector{2}) = a[1]*b[2]-a[2]*b[1]
cross(i,a::SVector{3},b::SVector{3}) = (j=i%3+1;k=(i+1)%3+1;a[j]*b[k]-a[k]*b[j])

"""
    Qcriterion(I::CartesianIndex,u)
    
    Compute the Q-criterion for a 2D/3D velocity field
"""
function Qcriterion(I,u)
    J = ∇u(I,u)
    S,Ω = (J+J')/2,(J-J')/2
    0.5*(√(tr(Ω*Ω'))^2-√(tr(S*S'))^2)
end
∇u(I::CartesianIndex{2},u) = @SMatrix [WaterLily.∂(i,j,I,u) for i ∈ 1:2, j ∈ 1:2]
∇u(I::CartesianIndex{3},u) = @SMatrix [WaterLily.∂(i,j,I,u) for i ∈ 1:3, j ∈ 1:3]

"""
    ∫2Qϕ!(::ForcePartitionMethod,::Flow,tᵢ,i,recompute=true)
    
    Compute the vorticity influence on the ith component of the force 
    on the body at a time tᵢ.
    - FPM: ForcePartitionMethod
    - a: Flow
    - tᵢ: time, the default is the current time
    - i: component of the force
    - recompute: if true, the potential is recomputed (needed for moving geometries)
"""
function ∫2QϕdV!(FPM::ForcePartitionMethod,a::Flow,tᵢ=sum(a.Δt);
                 axis=nothing,recompute=true,T=promote_type(Float64,eltype(a.p)))
    FPM.σ .= 0
    # get potential
    recompute && potential!(FPM,FPM.body;tᵢ=tᵢ,axis=axis)
    # compute the influence of the Q field
    @inside FPM.σ[I] = FPM.ϕ[I]*Qcriterion(I,a.u)
    # return the integral over the doman
    2sum(T,FPM.σ)
end

"""
    ∮UϕdS!(::ForcePartitionMethod,::Flow,tᵢ,i,recompute=true)

    Compute the kinematic influence on the `ith` component of the force.
    - FPM: ForcePartitionMethod
    - a: Flow
    - tᵢ: time, the default is the current time
    - i: component of the force
    - recompute: if true, the potential is recomputed (needed for moving geometries)
"""
function ∮UϕdS!(FPM::ForcePartitionMethod,a::Flow,tᵢ=sum(a.Δt);
                axis=nothing,recompute=true,T=promote_type(Float64,eltype(a.p)))
    FPM.σ .= 0
    # get potential
    recompute && potential!(FPM,FPM.body;tᵢ=tᵢ,axis=axis)
    # compute the influence of kinematics
    @inside FPM.σ[I] = FPM.ϕ[I]*dUdtnds(I,FPM.body,tᵢ)
    # return the integral over the body
    sum(T,FPM.σ)
end
using ForwardDiff: derivative
@inline function dUdtnds(I,body,tᵢ) #TODO check order of operation
    d,nᵢ,_ = measure(body,loc(0,I),tᵢ)
    aᵢ = derivative(tᵢ->derivative(tᵢ->body.map(loc(0,I),tᵢ),tᵢ),tᵢ)
    sum(nᵢ.*aᵢ)*WaterLily.kern(clamp(d,-1,1))
end
"""
    ∮ReωdS!(::ForcePartitionMethod,::Flow,tᵢ,i,recompute=true,type=:force)

    Compute the viscous influence on the `ith` component of the force. By specifying
    the `type` as `:moment` the influence on the moment is computed.
"""
function ∮ReωdS!(FPM::ForcePartitionMethod{A,T},a::Flow{D},tᵢ=sum(a.Δt);
                 axis=nothing,recompute=true) where {A,T,D}
    FPM.σ .= 0
    # get potential
    recompute && potential!(FPM,FPM.body;tᵢ=tᵢ,axis=axis)
    isnothing(axis) && (axis = A)
    e₁ = zeros(D); e₁[(axis-1)%3+1] = 1; e₁ = SVector{D}(e₁)
    # compute the vorticity × normal 
    @WaterLily.loop FPM.σ[I] = a.ν*dot(ωxn(I,a.u,FPM.body,tᵢ),∇ϕ(I,FPM.ϕ).-e₁) over I ∈ inside(FPM.σ)
    # return the integral over the body
    sum(promote_type(Float64,T),FPM.σ)
end
∇ϕ(I::CartesianIndex{2},ϕ) = @SVector[WaterLily.∂(i,I,ϕ) for i ∈ 1:2]
∇ϕ(I::CartesianIndex{3},ϕ) = @SVector[WaterLily.∂(i,I,ϕ) for i ∈ 1:3]
function ωxn(I::CartesianIndex{2},u,body,tᵢ)
    d,n,_ = measure(body,loc(0,I),tᵢ)
    ω = WaterLily.∂(2,1,I,u)-WaterLily.∂(1,2,I,u) # at cell center
    SA[-ω*n[2],ω*n[1]]*WaterLily.kern(clamp(d,-1,1))
end
import WaterLily: fSV,permute
cross(a::SVector{3},b::SVector{3}) = fSV(i->permute((j,k)->a[j]*b[k],i),3)
function ωxn(I::CartesianIndex{3},u,body,tᵢ)
    d,n,_ = measure(body,loc(0,I),tᵢ)
    ω = WaterLily.ω(I,u) # at cell center
    cross(ω,n)*WaterLily.kern(clamp(d,-1,1))
end