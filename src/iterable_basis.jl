abstract type Iterable <: ITensors.InterpolationType

struct IterableBasisInterpolation{T,N,TCoefs,IT<:DimSpec{Iterable},K<:Tuple{Vararg{AbstractVector}}} <: AbstractInterpolation{T,N,IT}
    knots::K
    coefs::TCoefs
    it::IT
end

"""
     Contsruct from a grid from. 
"""
function IterableBasisInterpolation(knots::Tuple{Vararg{...}}, A::AbstractArray, it::IT) where IT<:DimSpec{Iterable}
    interpolate(tweight(A), tcoef(A), knots, A, it)
end
function interpolate(::Type{Tweights}, ::Type{Tcoefs}, ...)
    interpolate(TWeights, knots, copy(A), it)
end
function interpolate!(::Type{Tweights}, ...)
    IterableBasisInterpolation(TWeights, knots, A, it)
end

interpolate!(values,it::IT) where IT<:DimSpec{Chebyshev}
@inline function (itp::IterableBasisInterpolation{T,N})(x::Vararg{Number,N}) where {T,N}
    wis = weightedindexes((value_weights,), itpinfo(itp)..., x)
    InterpGetindex(itp)[wis...]
end

@forward IterableBasisInterpolation.coefs Base.size, Base.axes
itpflag(A::IterableBasisInterpolation) = A.it
Base.parent(A::IterableBasisInterpolation) = A.coefs
coefficients(A::IterableBasisInterpolation) = A.coefs
