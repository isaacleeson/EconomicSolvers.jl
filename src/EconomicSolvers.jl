module EconomicSolvers

using CUDA, GPUArrays, Interpolations, ITensors
import Base: *, +, map!, copyto!

include("rowenhorst.jl")
include("map.jl")
include("transition.jl")
include("univariate_solvers.jl")

end
