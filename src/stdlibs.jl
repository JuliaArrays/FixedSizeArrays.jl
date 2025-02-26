using Random: Random

if isdefined(Base, :Memory)
    Random.rand!(rng::Random.AbstractRNG, A::FixedSizeArray{T,N,Memory{T}}, sp::Random.Sampler) where {N,T} =
        Random.rand!(rng, A.mem, sp)
end
