export norm_copy!

# map! with tuples as input and output arrays, and optional arguments
function Base.map!(g, dest::AnyGPUArray, input::NTuple{N,<:AnyGPUVector}, args...) where N
    map_kernel = (ctx, dest, f, input) -> begin
        cidx = GPUArrays.@cartesianidx dest
        val = f((@inbounds input[i][cidx[i]] for i in 1:N)...,args...)
        @inbounds dest[cidx] = val
        return nothing
    end
    gpu_call(map_kernel,dest,g,input)
    return dest
end
function Base.map!(g, dest::Tuple{Vararg{<:AnyGPUArray}}, input::NTuple{N,<:AnyGPUVector}, args...) where N
    map_kernel = (ctx, dest, f, input, dest_tail) -> begin
        cidx = GPUArrays.@cartesianidx dest
        val = f((@inbounds input[i][cidx[i]] for i in 1:N)...,vargs...)
        @inbounds dest[1][cidx] = val[1]
        for i in 2:N
            @inbounds dest[i][cidx] = val[i]
        end
        return nothing
    end
    gpu_call(map_kernel,dest[1],g,input,Base.tail(dest))
    return dest
end
function Base.map!(f, dest::AbstractArray, input::Tuple{Vararg{AbstractVector}}, args...)
    for cindx in CartesianIndices(Iterators.product(input...))
        val = f((@inbounds input[i][cidx[i]] for i in 1:N)...,args...)
        @inbounds dest[cidx] = val
    end
end
function Base.map!(f, dest::Tuple{Vararg{AbstractArray}}, input::Tuple{Vararg{AbstractVector}}, args...)
    for cindx in CartesianIndices(Iterators.product(input...))
        val = f((@inbounds input[i][cidx[i]] for i in 1:N)...,args...)
        @inbounds broadcast!(setindex!(_1[cidx],_2), dest, val)
    end
end

# norm copy
function norm_copy!(itr1,itr2,p::Real=2.0)
    err = mapreduce((x,y)->(x-y)^p,+,itr1,itr2)^(1/p)
    copyto!(itr1,itr2)
    return err
end

