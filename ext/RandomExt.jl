module RandomExt

using Random: Random
using FixedSizeArrays: FixedSizeArray

Random.rand!(rng::Random.AbstractRNG, A::FixedSizeArray{T,N}, sp::Random.Sampler) where {N,T} =
    Random.rand!(rng, parent(A), sp)

# Methods needed only to resolve ambiguities, we don't care too much about them
Random.rand!(r::Random.MersenneTwister, A::FixedSizeArray{Float64,N}, I::Random.SamplerTrivial{<:Random.FloatInterval{Float64}}) where {N} =
    Random.rand!(r, parent(A), I)
if VERSION < v"1.11"
    Random.rand!(::Random._GLOBAL_RNG, A::FixedSizeArray{Float64,N}, I::Random.SamplerTrivial{<:Random.FloatInterval{Float64}}) where {N} =
        Random.rand!(Random.default_rng(), parent(A), I)
end

end # module RandomExt
