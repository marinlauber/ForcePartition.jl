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
@inline center(i,I,u) = @inbounds (u[I,i]+u[I+δ(i,I),i])*0.5 # convert face to cell center
function span_average!(avrg::AverageFlow,flow::Flow)
    ϵ = inv(size(inside(flow.p),3)) # sum over domain
    N,n = size_u(flow.u); avrg.v .= 0.0; avrg.τ .= 0.0 # reset
    for i ∈ 1:2 # u = U - ũ
        # spanwise average of velocity component `i` at cell face
        @loop avrg.v[Base.front(I),i] += ϵ*flow.u[I,i] over I in inside(flow.p)
        # spanwise fluctuating velocity over all cells, needs to be cell centered
        @loop flow.f[I,i] = center(i,I,flow.u) - center(i,Base.front(I),avrg.v) over I in inside(flow.p)
    end 
    for i ∈ 1:2, j ∈ 1:2 # τᵢⱼᴿ = < u' ⊗ u' > , needs to be at cell center as ∇⋅τ is then added to faces
        @loop avrg.τ[Base.front(I),i,j] += ϵ*flow.f[I,i].*flow.f[I,j] over I in inside(flow.p)
    end
end

# fluctuating velocity u' = u - U at cell center I
fluct(a,I,u,U) = @inbounds 0.5*((u[I+δ(a,I),a]+u[I,a])-(U[I+δ(a,I),a]+U[I,a]))
function mean!(avrg::AverageFlow, flow::Flow; stats_turb=true)
    dt = time(flow) - avrg.t[end]
    ϵ = dt / (dt + (avrg.t[end] - avrg.t[1]) + eps(eltype(flow.p)))
    @loop avrg.ϕ[I] = ϵ*flow.p[I] + (1.0 - ϵ)*avrg.ϕ[I] over I in CartesianIndices(flow.p)
    @loop avrg.v[Ii] = ϵ*flow.u[Ii] + (1.0 - ϵ)*avrg.v[Ii] over Ii in CartesianIndices(flow.u)
    if stats_turb
        for i in 1:ndims(flow.p), j in 1:ndims(flow.p) # τᵢⱼᴿ = < u' ⊗ u' > , needs to be at cell center as ∇⋅τ is then added to faces
            @loop avrg.τ[I,i,j] = ϵ*(fluct(i,I,flow.u,avrg.v)*fluct(j,I,flow.u,avrg.v))+(1.0-ϵ)*avrg.τ[I,i,j] over I in inside(flow.p)
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

function spread!(src::Flow{2}, dest::Flow{3}; ϵ=0)
    @assert size(src.p)==Base.front(size(dest.p)) "a::Flow{2} must be the same size as b::Flow{3}[:,:,1,i] to spread"
    destp=dest.p; srcp=src.p; destu=dest.u; srcu=src.u # alias
    @loop destp[I] = srcp[Base.front(I)] over I in inside(destp)
    for i ∈ 1:2 # can only spread 2 components
        @loop destu[I,i] = srcu[Base.front(I),i]+ϵ*rand() over I in inside(destp)
    end
end
