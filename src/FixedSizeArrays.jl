module FixedSizeArrays

export FixedSizeArray, FixedSizeVector, FixedSizeMatrix

struct FixedSizeArray{T,N} <: DenseArray{T,N}
    mem::Memory{T}
    size::NTuple{N,Int}
end

const FixedSizeVector{T} = FixedSizeArray{T,1}
const FixedSizeMatrix{T} = FixedSizeArray{T,2}

function (self::Type{FixedSizeArray{T,N}})(::UndefInitializer, size::Vararg{Int,N}) where {T,N}
    return FixedSizeArray(Memory{T}(undef, prod(size)), size)
end

Base.IndexStyle(::Type{<:FixedSizeArray}) = IndexLinear()
Base.@propagate_inbounds Base.getindex(A::FixedSizeArray, i::Int) = A.mem[i]
Base.@propagate_inbounds Base.setindex!(A::FixedSizeArray, v, i::Int) = A.mem[i] = v

Base.size(a::FixedSizeArray) = getfield(a, :size)

end # module FixedSizeArrays
