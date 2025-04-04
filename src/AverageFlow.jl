using WaterLily
using StatsBase
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

# f = ∇⋅τ = ∂ᵢτᵢⱼ, for example f₁ = ∂₁τ₁₁ + ∂₂τ₂₁, because τ is at cell center inline components are
# simply the difference between neighboring cells, but the diagonal components are intepolated using the
# average of the two neighboring cells.
@fastmath SANS!(f::AbstractArray,τ::AbstractArray) = (f.=0; for (i,j) ∈ Iterators.product(1:2,1:2)
    @loop f[I,i] = (i==j ? τ[I,i,j]-τ[I-δ(i,I),i,j] :
                    @inbounds(τ[I+δ(j,I),i,j]+τ[I-δ(i,I)+δ(j,I),i,j]
                             -τ[I-δ(j,I),i,j]-τ[I-δ(i,I)-δ(j,I),i,j])/4)  over I in inside_u(f)
end)
"""
    spread!(src::Flow{2}, dest::Flow{3}; ϵ=0)

    Spread the 2D flow field to a 3D flow field. The 2D flow field is spread to the 3rd dimension of the 3D flow field.
    src - 2D flow field
    dest - 3D flow field
    ϵ - random noise added to the spread
"""
function spread!(src::Flow{2}, dest::Flow{3}; ϵ=0)
    @assert size(src.p)==Base.front(size(dest.p)) "a::Flow{2} must be the same size as b::Flow{3}[:,:,1,i] to spread"
    destp=dest.p; srcp=src.p; destu=dest.u; srcu=src.u # alias
    @loop destp[I] = srcp[Base.front(I)] over I in inside(destp)
    for i ∈ 1:2 # can only spread 2 components
        @loop destu[I,i] = srcu[Base.front(I),i]+ϵ*rand() over I in inside(destp)
    end
end

"""
    profile!(sim,profiles;loc=[1,2])

    Measure velocity profiles at certain locations.
    sim - simulation object
    profiles - array of profiles
    loc - locations to measure the profiles in lengthscale units

"""
function profile!(sim,profiles;loc=[1,2]) # measure velocity profiles at certain locations
    l = []
    for x ∈ loc.*sim.L
        push!(l,azimuthal_avrg(sim.flow.u[Int(x),:,:,1]))
    end
    push!(profiles,l)
end
Base.hypot(I::CartesianIndex) = √sum(abs2,I.I)
"""
    azimuthal_avrg(data; center=nothing, binsize=1.0)

    Calculate the azimuthally averaged radial profile.
    image - The 2D image
    center - The [x,y] pixel coordinates used as the center. The default is 
             None, which then uses the center of the image (including 
             fractional pixels).
    binsize - size of the averaging bin.  Can lead to strange results if
        non-binsize factors are used to specify the center and the binsize is
        too large
    """
function azimuthal_avrg(data; center=nothing, binsize=1.0)

    # Calculate the indices from the image
    CIs = CartesianIndices(data)
    isnothing(center) && (center = (maximum(CIs)-minimum(CIs)).I.÷2 .+ 1)
    
    # radial distance from the center, make it a vector
    r = weights(hypot.(collect(CIs .- CartesianIndex(center)))).values
    
    # the 'bins' as initially defined are lower/upper bounds for each bin
    # so that values will be in [lower,upper)  
    nbins = Int(round(maximum(r) / binsize))
    maxbin = nbins * binsize
    bins = range(0,maxbin;length=nbins)
    # but we're probably more interested in the bin centers than their left or right sides...
    bin_centers = (bins[1:end-1].+bins[2:end])/2.0
    r_weights = fit(Histogram, r, bins, closed=:left).weights

    # compute the azimuthal average
    radial_prof = fit(Histogram, r, weights(data), bins, closed=:left).weights ./ r_weights
    
    return bin_centers, radial_prof
end