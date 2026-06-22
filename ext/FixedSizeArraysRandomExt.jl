module FixedSizeArraysRandomExt

using Random: Random
using FixedSizeArrays: FixedSizeArray

function Random.rand!(rng::Random.AbstractRNG, A::FixedSizeArray{T,N}, sp::Random.Sampler) where {N,T}
    Random.rand!(rng, parent(A), sp)
    return A
end

# Methods needed only to resolve ambiguities, we don't care too much about them
function Random.rand!(r::Random.MersenneTwister, A::FixedSizeArray{Float64,N}, I::Random.SamplerTrivial{<:Random.FloatInterval{Float64}}) where {N}
    Random.rand!(r, parent(A), I)
    return A
end

if VERSION < v"1.11"
    function Random.rand!(::Random._GLOBAL_RNG, A::FixedSizeArray{Float64,N}, I::Random.SamplerTrivial{<:Random.FloatInterval{Float64}}) where {N}
        Random.rand!(Random.default_rng(), parent(A), I)
        return A
    end
end

# UnsafeRefArray-backed: delegate rand! to the underlying Memory so that the
# specialized SIMD path in Random.XoshiroSimd is hit.
if @isdefined Memory
    using FixedSizeArrays: UnsafeRefArray

    function Random.rand!(rng::Random.AbstractRNG, A::UnsafeRefArray{T}, sp::Random.Sampler) where {T}
        Random.rand!(rng, parent(A), sp)
        return A
    end

    function Random.rand!(r::Random.MersenneTwister, A::UnsafeRefArray{Float64}, I::Random.SamplerTrivial{<:Random.FloatInterval{Float64}})
        Random.rand!(r, parent(A), I)
        return A
    end
end

end # module FixedSizeArraysRandomExt
