export transition!

function transition!(output,input,rule::Interpolations.InterpolationType)
    for cidx in CartesianIndices(A)
        idx_weights = Interpolations.weightedindexes((Interpolations.value_weights,),Interpolations.itpinfo(rule)...,Tuple(cidx))
        it = zip(Iterators.product(Interpolations.indextuple.(idx_weights)...), Iterators.product(Interpolations.weights.(idx_weights)...))
        for (idx, weight) in it
            @inbounds output[idx...] = prod(weight)*input[idx...]
        end
    end
end

