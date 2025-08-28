# impulse force due to vorticityusing WaterLily,StaticArrays,CUDA,BiotSavartBCs
struct VortexImpulse{D,T,Vf<:AbstractArray{T}}
    ω  :: Vf
    Ix :: Array{T}
    time :: Array{T}
    function VortexImpulse(sim;axis=1)
        ω = sim.flow.f # we point to the flow arrays
        # initial impulse
        @inside ω[I] = WaterLily.curl(3,I,sim.flow.u)*WaterLily.loc(2,I)[2]
        Ix = [sum(ωy[inside(ωy)])]
        time = [sim_time(sim)]
        # return struct
        new{ndims(ω)-1,typeof(ω),eltype(Ix)}(ω,Ix,time)
    end
end

import WaterLily: @loop
function vortex_force(Iω::VortexImpulse{2},sim;axis=1)
    j,sgn = axis%2+1,(-1)^(axis+1) # 1->y, 2->x
    @loop Iω.ω[I,1] = sgn*WaterLily.curl(3,I,sim.flow.u)*WaterLily.loc(j,I)[j]
    push!(Iω.time, sim_time(sim))
    push!(Iω.Ix, sum(ωy[inside(ω)]))
    return (Iω.Ix[end]-Iω.Ix[end-1])/(Iω.time[end]-Iω.time[end-1])
end
function vortex_force(Iω::VortexImpulse{3},sim;axis=1)
    j,sgn = axis%2+1,(-1)^(axis+1) # 1->y, 2->x
    @loop Iω.ω[Ii] = WaterLily.curl(Ii,sim.flow.u)*WaterLily.loc(2,I)[2] over I in inside_u(Iω.ω)
    push!(Iω.time,sim_time(sim))
    push!(Iω.Ix, sum(ωy[inside(ωy)]))
    return (Iω.Ix[end]-Iω.Ix[end-1])/(Iω.time[end]-Iω.time[end-1])
end