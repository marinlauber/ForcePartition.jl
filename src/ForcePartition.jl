module ForcePartition

include("ForcePartitionMethod.jl")
export ForcePartitionMethod,potential!,∫2QϕdV!,∮UϕdS!,∮ReωdS!

include("AverageFlow.jl")
export AverageFlow,SpanAverage,MeanFlow,span_average!,mean!,write!,spread!,SANS!

include("VorticityForce.jl")
export VortexImpulse

end
