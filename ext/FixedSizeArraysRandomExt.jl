module FixedSizeArraysRandomExt

using Random: Random
using FixedSizeArrays: FixedSizeArrayDims

function Random.rand!(rng::Random.AbstractRNG, A::FixedSizeArrayDims{T}, sp::Random.Sampler) where {T}
    Random.rand!(rng, parent(A), sp)
    return A
end

# Methods needed only to resolve ambiguities, we don't care too much about them
function Random.rand!(r::Random.MersenneTwister, A::FixedSizeArrayDims{Float64}, I::Random.SamplerTrivial{<:Random.FloatInterval{Float64}})
    Random.rand!(r, parent(A), I)
    return A
end

if VERSION < v"1.11"
    function Random.rand!(::Random._GLOBAL_RNG, A::FixedSizeArrayDims{Float64}, I::Random.SamplerTrivial{<:Random.FloatInterval{Float64}})
        Random.rand!(Random.default_rng(), parent(A), I)
        return A
    end
end

end # module FixedSizeArraysRandomExt
