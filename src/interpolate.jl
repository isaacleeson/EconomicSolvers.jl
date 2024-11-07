function interpolate!(ranges::Tuple{Vararg{UnitRange}}, storage::AbstractArray, interpolation_type::Interpolations.DimSpec{Gridded})
    return interpolate!(storage,interpolation_type)
end
function interpolate!(range::Tuple{Vararg{AbstractRange}}, storage::AbstractArray, interpolation_type::Interpolations.DimSpec{BSpline})
    return scale(interpolate!(storage,interpolation_type),ranges)
end
function copyto!(itp::Interpolations.BSplineInterpolation, values::AbstractArray)
    Interpolations.coefficients(itp) .= values
    interpolate!(Interpolations.knots(itp), Interpolations.coefficients(itp), Interpolations.itpflag(itp))
end
function copyto!(itp::Interpolations.ScaledInterpolation, values::AbstractArray)
    Interpolations.coefficients(itp) .= values
    return scale(interpolate!(Interpolations.knots(itp), Interpolations.coefficients(itp),itp.ranges))
end
function copyto!(itp::Interpolations.GriddedInterpolation, values::AbstractArray)
    Interpolations.coefficients(itp) .= values
    interpolate!(itp.knots, Interpolations.coefficients(itp), Interpolations.itpflag(itp))
end
