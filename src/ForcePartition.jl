module ForcePartition

using WaterLily: AutoBody
using ParametricBodies: ParametricBody, HashedLocator

# needed for grabbing the bodies
Base.copy(b::AutoBody) = AutoBody(b.sdf,b.map)
Base.copy(b::ParametricBody) = ParametricBody(b.surf,copy(b.locate),b.map,b.scale)
Base.copy(l::HashedLocator) = HashedLocator(l.refine,l.lims,copy(l.hash),l.lower,l.step)

include("ForcePartitionMethod.jl")
export ForcePartitionMethod, potential!, ∫2Qϕ!

end
