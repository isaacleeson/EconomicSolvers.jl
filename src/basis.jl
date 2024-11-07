struct Basis{F, N, AR, IT, AX} <: ITensors.AbstractInterpolation{F, N, IT}
end
struct Chebyshev{D} <: ITensors.InterpolationType
    degree:D
end
