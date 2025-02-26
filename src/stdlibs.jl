using Random: Random

if isdefined(Base, :Memory)
    Random.rand!(rng::Random.AbstractRNG, A::FixedSizeArray{T,N,Memory{T}}, sp::Random.Sampler) where {N,T} =
        Random.rand!(rng, A.mem, sp)

    # Method needed only to resolve ambiguities, we don't care too much about it.
    Random.rand!(r::Random.MersenneTwister, A::FixedSizeArrayDefault{Float64,N}, I::Random.SamplerTrivial{<:Random.FloatInterval{Float64}}) where {N} =
        Random.rand!(r, A.mem, I)
end
