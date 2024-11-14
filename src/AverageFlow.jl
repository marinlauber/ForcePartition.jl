using WaterLily
using JLD2

"""
     AverageFlow{T, Sf<:AbstractArray{T}, Vf<:AbstractArray{T}, Mf<:AbstractArray{T}}
Holds averages of velocity, , pressure, and Reynolds stresses.
"""
struct AverageFlow{T, Sf<:AbstractArray{T}, Vf<:AbstractArray{T}, Mf<:AbstractArray{T}}
    ϕ :: Sf # scalar, cell senter
    v :: Vf # vector, face
    τ :: Mf # tensor, cell center
    t :: Vector{T} # time
end
function SpanAverage(flow::Flow{D,T}; N=size(flow.p), t_init=0.0) where {D,T}
    mem = typeof(flow.u).name.wrapper
    ϕ = zeros(T, Base.front(N)) |> mem
    v = zeros(T, Base.front(N)...,D-1) |> mem
    τ = zeros(T, Base.front(N)...,D-1,D-1) |> mem
    AverageFlow(ϕ,v,τ,T[t_init])
end
function MeanFlow(flow::Flow{D,T}; N=size(flow.p), t_init=0.0) where {D,T}
    mem = typeof(flow.u).name.wrapper
    ϕ = zeros(T, N) |> mem
    v = zeros(T, N...,D) |> mem
    τ = zeros(T, N...,D,D) |> mem
    AverageFlow(ϕ,v,τ,T[t_init])
end

import WaterLily: Flow,@loop,time,size_u,inside_u,inside
# spanwise average of velocity component `a` at cell center
∫dz(a,I,u,N) = @inbounds sum(@views(u[I,2:N[3]-1,a]))
center(i,I,u) = @inbounds (u[I,i]+u[I+δ(i,I),i])*0.5 # convert face to cell center
function span_average!(avrg::AverageFlow,flow::Flow)
    ϵ = inv(size(flow.p,3)-2) # sum over domain
    N,n = size_u(flow.u)
    for i ∈ 1:2 # u = U - ũ
        # spanwise average of velocity component `i` at cell face
        @loop avrg.v[I,i] = ϵ*∫dz(i,I,flow.u,N) over I in inside(avrg.ϕ)
        # spanwise fluctuating velocity over all cells, needs to be cell centered
        @loop flow.f[I,i] = center(i,I,flow.u) - center(i,Base.front(I),avrg.v) over I in inside(flow.p)
    end 
    for i ∈ 1:2, j ∈ 1:2 # τᵢⱼᴿ = < u' ⊗ u' > , needs to be at cell center as ∇⋅τ is then added to faces
        @loop avrg.τ[I,i,j] = ϵ*sum(flow.f[I,2:N[3]-1,i].*flow.f[I,2:N[3]-1,j]) over I in CartesianIndices(avrg.ϕ)
    end
end

# fluctuating velocity u' = u - U at cell center I
fluct(a,I,u,U) = @inbounds 0.5*((u[I+δ(a,I),a]+u[I,a])-(U[I+δ(a,I),a]+U[I,a]))
function mean!(avrg::AverageFlow, flow::Flow; stats_turb=true)
    dt = time(flow) - avrg.t[end]
    ε = dt / (dt + (avrg.t[end] - avrg.t[1]) + eps(eltype(flow.p)))
    @loop avrg.ϕ[I] = ε*flow.p[I] + (1.0 - ε)*avrg.ϕ[I] over I in CartesianIndices(flow.p)
    @loop avrg.v[Ii] = ε*flow.u[Ii] + (1.0 - ε)*avrg.v[Ii] over Ii in CartesianIndices(flow.u)
    if stats_turb
        for i in 1:ndims(flow.p), j in 1:ndims(flow.p) # τᵢⱼᴿ = < u' ⊗ u' > , needs to be at cell center as ∇⋅τ is then added to faces
            @loop avrg.τ[I,i,j] = ϵ*(fluct(i,I,flow.u,avrg.v)*fluct(j,I,flow.u,avrg.v))+(1.0-ε)*avrg.τ[I,i,j] over I in inside(flow.p)
        end
    end
    push!(avrg.t, avrg.t[end] + dt)
end

WaterLily.write!(fname, avrgflow::AverageFlow; dir="data/") = jldsave(
    dir*fname*".jld2";
    ϕ=Array(avrgflow.ϕ),
    f=Array(avrgflow.f),
    τ=Array(avrgflow.τ),
    t=avrgflow.t
)

# f = ∇⋅τ = ∂ᵢτᵢⱼ
SANS!(f::AbstractArray,τ::AbstractArray) = (f.=0; for (i,j) ∈ Iterators.product(1:2,1:2)
    @loop f[I,i] = τ[I,i,j]-τ[I-δ(i,I),i,j] over I in inside_u(f)
end)

import WaterLily: scale_u!,BCTuple,conv_diff!,accelerate!,BDIM!,BC!,project!,CFL
@fastmath function mom_step_force!(a::Flow{N},b::AbstractPoisson,∇τ) where N
    a.u⁰ .= a.u; scale_u!(a,0)
    # predictor u → u'
    U = BCTuple(a.U,@view(a.Δt[1:end-1]),N)
    conv_diff!(a.f,a.u⁰,a.σ,ν=a.ν,perdir=a.perdir)
    accelerate!(a.f,@view(a.Δt[1:end-1]),a.g,a.U)
    force!(a,∇τ) # new forcing
    BDIM!(a); BC!(a.u,U,a.exitBC,a.perdir)
    a.exitBC && exitBC!(a.u,a.u⁰,U,a.Δt[end]) # convective exit
    project!(a,b); BC!(a.u,U,a.exitBC,a.perdir)
    # corrector u → u¹
    U = BCTuple(a.U,a.Δt,N)
    conv_diff!(a.f,a.u,a.σ,ν=a.ν,perdir=a.perdir)
    accelerate!(a.f,a.Δt,a.g,a.U)
    force!(a,∇τ) # new forcing
    BDIM!(a); scale_u!(a,0.5); BC!(a.u,U,a.exitBC,a.perdir)
    project!(a,b,0.5); BC!(a.u,U,a.exitBC,a.perdir)
    push!(a.Δt,CFL(a))
end
force!(a::Flow,Ⅎ::AbstractArray) = @loop a.f[Ii] += Ⅎ[Ii] over Ii in CartesianIndices(a.f)

function spread!(src::Flow{2}, dest::Flow{3}; ϵ=0)
    @assert size(src.p)==Base.front(size(dest.p)) "a::Flow{2} must be the same size as b::Flow{3}[:,:,1,i] to spread"
    destp=dest.p; srcp=src.p; destu=dest.u; srcu=src.u # alias
    @loop destp[I] = srcp[Base.front(I)] over I in inside(destp)
    for i ∈ 1:2 # can only spread 2 components
        @loop destu[I,i] = srcu[Base.front(I),i]+ϵ*rand() over I in inside(destp)
    end
end
