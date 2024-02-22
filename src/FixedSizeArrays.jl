module FixedSizeArrays

export FixedSizeArray, FixedSizeVector, FixedSizeMatrix

struct FixedSizeArray{T,N} <: DenseArray{T,N}
    mem::Memory{T}
    const size::NTuple{N,Int}
end

const FixedSizeVector{T} = FixedSizeArray{T,1}
const FixedSizeMatrix{T} = FixedSizeArray{T,2}

function (self::Type{FixedSizeArray{T,N}})(::UndefInitializer, size::Vararg{Int,N}) where {T,N}
    return FixedSizeArray(Memory{T}(undef, prod(size)), size)
end

Base.@propagate_inbounds function Base.setindex!(A::FixedSizeArray{T}, x, i::Int) where {T}
    getfield(A, :mem)[i] = x
    return A
end
Base.@inline function Base.setindex!(A::FixedSizeArray{T}, x, i1::Int, i2::Int, I::Int...) where {T}
    @boundscheck checkbounds(A, i1, i2, I...) # generally _to_linear_index requires bounds checking
    getfield(A, :mem)[Base._to_linear_index(A, i1, i2, I...)] = x
    return A
end
Base.@propagate_inbounds function Base.getindex(A::FixedSizeArray, i::Int)
    getfield(A, :mem)[i]
end
function Base.getindex(A::FixedSizeArray, i1::Int, i2::Int, I::Int...)
    @inline
    @boundscheck checkbounds(A, i1, i2, I...) # generally _to_linear_index requires bounds checking
    return @inbounds A[Base._to_linear_index(A, i1, i2, I...)]
end

Base.size(a::FixedSizeArray) = getfield(a, :size)

end # module FixedSizeArrays
