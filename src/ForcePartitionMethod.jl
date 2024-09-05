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
struct ForcePartitionMethod{T,Sf<:AbstractArray{T},Vf<:AbstractArray{T}}
    ϕ :: Sf
    σ :: Sf
    f :: Vf
    pois :: AbstractPoisson
    body :: AbstractBody
    function ForcePartitionMethod(sim::Simulation)
        ϕ,σ,f = copy(sim.flow.p),sim.flow.σ,sim.flow.f # σ and f point to the flow fields
        new{eltype(sim.flow.u),typeof(ϕ),typeof(f)}(ϕ,σ,f,sim.pois,sim.body) # we point to the poisson solver of the flow
    end
end

test(a,::Val{:force}) = println("force")
test(a,::Val{:moment}) = println("moment")
func(a,f::Symbol=:force) = test(a,Val(f))
"""
    potential!(FPM::ForcePartitionMethod;axis,tᵢ,f::Symbol=:force)

    Compute the the influence potential of the body at time tᵢ for either 
    the `:force` or the `:moment`.
"""
function potential!(FPM::ForcePartitionMethod{T};x₀=0,tᵢ=0,axis=1,type::Symbol=:force) where T
    # first we need to save the pressure field
    @inside FPM.σ[I] = FPM.pois.x[I]
    # generate source term
    @inside FPM.pois.z[I] = source(FPM.body,loc(0,I,T),x₀,axis,tᵢ,Val(type))
    # solver for potential
    solver!(FPM.pois) # keep the tol the same as the pressure
    # copy to the FPM
    @inside FPM.ϕ[I] = FPM.pois.x[I]
    @inside FPM.pois.x[I] = FPM.σ[I] # put back pressure field
end

"""
    force(body,x₀,i,tᵢ,::Val{:force})
    
    Compute the the source term for the force influence potential
    of the body at time tᵢ
"""
function source(body,x,x₀,i,tᵢ,::Val{:force})
    dᵢ,nᵢ,_ = measure(body,x,tᵢ)
    WaterLily.kern(clamp(dᵢ,-1,1))*nᵢ[i] # i-component of the normal
end
"""
    moment(body,x,x₀,i,tᵢ,::Val{:moment})

    Compute the the source term for the moment influence potential
    of the body at time tᵢ
"""
function source(body,x,x₀,i,tᵢ,::Val{:moment}) # the axis here does not matter
    dᵢ,nᵢ,_ = measure(body,x,tᵢ)
    WaterLily.kern(clamp(dᵢ,-1,1))*cross((x.-x₀),nᵢ) # the normal moment
end
cross(a::SVector{2},b::SVector{2}) = a[1]*b[2]-a[2]*b[1]

"""
    Qcriterion(I::CartesianIndex,u)
    
    Compute the Q-criterionfor 2D and 3D velocity fields
"""
function Qcriterion(I::CartesianIndex{2},u)
    J = @SMatrix [WaterLily.∂(i,j,I,u) for i ∈ 1:2, j ∈ 1:2]
    S,Ω = (J+J')/2,(J-J')/2
    ## -0.5*sum(eigvals(S^2+Ω^2)) # this is also possible, but 2x slower
    0.5*(√(tr(Ω*Ω'))^2-√(tr(S*S'))^2)
end
function Qcriterion(I::CartesianIndex{3},u)
    J = @SMatrix [WaterLily.∂(i,j,I,u) for i ∈ 1:3, j ∈ 1:3]
    S,Ω = (J+J')/2,(J-J')/2
    ## -0.5*sum(eigvals(S^2+Ω^2)) # this is also possible, but 2x slower
    0.5*(√(tr(Ω*Ω'))^2-√(tr(S*S'))^2)
end

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
function ∫2Qϕ!(FPM::ForcePartitionMethod,a::Flow,tᵢ=WaterLily.time(a);axis=1,recompute=true,type::Symbol=:force)
    # get potential
    recompute && potential!(FPM;type=type,tᵢ=tᵢ,axis=axis)
    # compute the influence of the Q field
    @inside FPM.σ[I] = FPM.ϕ[I]*Qcriterion(I,a.u)
    # return the integral over the doman
    2sum(@inbounds(FPM.σ[inside(FPM.σ)]))
end

function ∮UϕdS!(FPM::ForcePartitionMethod,a::Flow,tᵢ=WaterLily.time(a);axis=1,recompute=true,type::Symbol=:force)
    # get potential
    recompute && potential!(FPM;type=type,tᵢ=tᵢ,axis=axis)
    # compute the influence of kinematics
    @inside FPM.σ[I] = FPM.ϕ[I]*dUdtnds(I,FPM.body,tᵢ)
    # return the integral over the body
    sum(@inbounds(FPM.σ[inside(FPM.σ)]))
end
using ForwardDiff: derivative
@inline function dUdtnds(I,body,tᵢ)
    d,nᵢ,_ = measure(body,loc(0,I),tᵢ)
    aᵢ = derivative(tᵢ->derivative(tᵢ->body.map(loc(0,I),tᵢ),tᵢ),tᵢ)
    sum(nᵢ.*aᵢ)*WaterLily.kern(clamp(d,-1,1))
end

function ∮ReωdS!(FPM::ForcePartitionMethod{T},a::Flow{D},tᵢ=WaterLily.time(a);
                 axis=1,recompute=true,type::Symbol=:force) where {T,D}
    # get potential
    recompute && potential!(FPM;type=type,tᵢ=tᵢ,axis=axis)
    e₁ = zeros(D); e₁[axis] = 1; e₁ = SVector{D}(e₁)
    # compute the vorticity × normal 
    @WaterLily.loop FPM.σ[I] = dot(ωxn(I,a.u,FPM.body,tᵢ),∇ϕ(I,FPM.ϕ).-e₁) over I ∈ inside(FPM.σ)
    # return the integral over the body
    a.ν*sum(@inbounds(FPM.σ[inside(FPM.σ)]))
end
∇ϕ(I::CartesianIndex{2},ϕ) = @SVector[WaterLily.∂(i,I,ϕ) for i ∈ 1:2]
∇ϕ(I::CartesianIndex{3},ϕ) = @SVector[WaterLily.∂(i,I,ϕ) for i ∈ 1:3]
# @TODO vorticity at the cell center
function ωxn(I::CartesianIndex{2},u,body,tᵢ)
    d,n,_ = measure(body,loc(0,I),tᵢ)
    ω = WaterLily.curl(3,I,u) # this one is NOT at cell center
    SA[-ω*n[2],ω*n[1]]*WaterLily.kern(clamp(d,-1,1))
end
function ωxn(I,u,body,tᵢ)
    d,n,_ = measure(body,loc(0,I),tᵢ)
    ω = WaterLily.ω(I,u) # at cell center
    cross(ω,n)*WaterLily.kern(clamp(d,-1,1))
end