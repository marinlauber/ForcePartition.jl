using WaterLily
using StaticArrays
using LinearAlgebra: tr

struct ForcePartitionMethod{T,Sf<:AbstractArray{T},Vf<:AbstractArray{T}}
    ϕ :: Sf
    σ :: Sf
    f :: Vf
    pois::AbstractPoisson
    body::AbstractBody
    function ForcePartitionMethod(sim::Simulation)
        ϕ,σ,f = copy(sim.flow.p),sim.flow.σ,sim.flow.f # σ and f point to the flow fields
        new{eltype(sim.flow.u),typeof(ϕ),typeof(f)}(ϕ,σ,f,sim.pois,sim.body) # we point to the poisson solver of the flow
    end
end

"""
    Compute the the influence potential of the body at time tᵢ 
"""
# #@TODO check that updating the body updates the poisson coefficients of the FPM
function potential!(FPM::ForcePartitionMethod;x₀=SA[0.,0.],f::Symbol=:force,tᵢ=0,axis=1,tol=1e-5,itmx=1e4)
    # first we need to save the pressure field
    @inside FPM.σ[I] = FPM.pois.x[I]
    # generate source term
    apply!(x->eval(f)(FPM.body,x,x₀,axis,tᵢ),FPM.pois.z)
    # solver for potential
    solver!(FPM.pois;itmx=itmx,tol=tol)
    # copy to the FPM
    @inside FPM.ϕ[I] = FPM.pois.x[I]
    @inside FPM.pois.x[I] = FPM.σ[I] # put back pressure field
end

# the boundary value is the surface normal
function force(body,x,x₀,i,tᵢ)
    dᵢ,nᵢ,_ = measure(body,x,tᵢ)
    WaterLily.kern(clamp(dᵢ,-1,1))*nᵢ[i] # i-component of the normal
end
cross(a::SVector{2},b::SVector{2}) = a[1]*b[2]-a[2]*b[1]
# the boundary value is the surface normal
function moment(body,x,x₀,i,tᵢ) # the axis here does not matter
    dᵢ,nᵢ,_ = measure(body,x,tᵢ)
    WaterLily.kern(clamp(dᵢ,-1,1))*cross((x-x₀),nᵢ) # the normal moment
end

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

function ∫2Qϕ!(FPM::ForcePartitionMethod,a::Flow,tᵢ=WaterLily.time(a);axis=1,recompute=true,f::Symbol=:force)
    # get potential
    recompute && potential!(FPM;f=f,tᵢ=tᵢ,axis=axis)
    # compute the influence of the Q field
    @WaterLily.loop FPM.σ[I] = FPM.ϕ[I]*Qcriterion(I,a.u) over I ∈ inside(FPM.σ)
    # return the integral over the doman
    2sum(@inbounds(FPM.σ[inside(FPM.σ)]))
end

function ∮Ubϕ!(FPM::ForcePartitionMethod,a::Flow,tᵢ=WaterLily.time(a);axis=1,recompute=true,f::Symbol=:force)
    # get potential
    recompute && potential!(FPM;f=f,tᵢ=tᵢ,axis=axis)
    # compute the influence of kinematics
    @WaterLily.loop FPM.σ[I] = FPM.ϕ[I]*Unds(I,FPM.body,tᵢ) over I ∈ inside(FPM.σ)
    # return the integral over the doman
    sum(@inbounds(FPM.σ[inside(FPM.σ)]))
end
@inline function Unds(I,body,tᵢ)
    d,nᵢ,vᵢ = measure(body,loc(0,I),tᵢ)
    sum(nᵢ.*vᵢ)*WaterLily.kern(clamp(d,-1,1))
end

function ∮ReωdS!(FPM::ForcePartitionMethod,a::Flow,tᵢ=WaterLily.time(a);axis=1,recompute=true,f::Symbol=:force)
    # get potential
    recompute && potential!(FPM;f=f,tᵢ=tᵢ,axis=axis)
    # compute the vorticity × normal 
    for i ∈ 1:2
        @WaterLily.loop FPM.f[I,i] = ωxn(i,I,a.u,body,tᵢ)*WaterLily.∂(i,I,FPM.ϕ) over I ∈ inside(FPM.σ)
    end
    # compute the influence of the Q field
    @WaterLily.loop FPM.σ[I] = FPM.ϕ[I]*Unds(I,FPM.body,tᵢ) over I ∈ inside(FPM.σ)

end

function ωxn(i,I,u,body,tᵢ)
    d,n,_ = measure(body,loc(i,I),tᵢ)
    ω = WaterLily.curl(i,I,u)
    -sign(2-i)*ω*n[i%2+1]*WaterLily.kern(clamp(d,-1,1))
end