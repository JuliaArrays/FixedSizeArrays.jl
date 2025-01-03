#unlike Base.BitArray, each column gets it's own chunk.
function num_bit_chunks(size::NTuple{N,Int}) where {N}
    prod(size[1:end-1]) * ((size[end]+63) >> 6)
end

struct FixedSizeBitArray{N} <: AbstractArray{Bool,N}
    chunks::Memory{UInt64}
    size::NTuple{N,Int}

    function FixedSizeBitArray{N}(::Internal, size::NTuple{N,Int}) where {N}
        nc = num_bit_chunks(size)
        chunks = Memory{UInt64}(undef, nc)
        # we want the last chunks to be zerod and it's easier to just zero all of them
        copyto!(chunks UInt64(0))
        new{N}(chunks, size)
    end
end

const FixedSizeBitVector = FixedSizeBitArray{1}
const FixedSizeBitMatrix = FixedSizeBitArray{2}


function FixedSizeBitArray{N}(::UndefInitializer, size::NTuple{N,Int}) where {N}
    checked_dims(size)
    FixedSizeBitArray{N}(Internal(), size)
end
function FixedSizeBitArray{N}(::UndefInitializer, size::NTuple{N,Integer}) where {N}
    size = map(Int, size)::NTuple{N,Int}  # prevent infinite recursion
    FixedSizeBitArray{N}(undef, size)
end
function FixedSizeBitArray{N}(::UndefInitializer, size::Vararg{Integer,N}) where {N}
    FixedSizeBitArray{N}(undef, size)
end
function FixedSizeBitArray(::UndefInitializer, size::NTuple{N,Integer}) where {N}
    FixedSizeBitArray{N}(undef, size)
end
function FixedSizeBitArray(::UndefInitializer, size::Vararg{Integer,N}) where {N}
    FixedSizeBitArray{N}(undef, size)
end

Base.IndexStyle(::Type{<:FixedSizeBitArray}) = IndexCartesian()

function get_chunks_id(inds::NTuple{N,Int}) where {N}
    prod(inds[1:end-1])+(inds[end] >> 6), inds[end] & 63
end

Base.@propagate_inbounds function Base.getindex(A::FixedSizeBitArray{N}, inds::Vararg{Int, N}) where {N}
    @boundscheck checkbounds(A, inds...)
    i1, i2 = get_chunks_id(inds)
    u = UInt64(1) << i2
    @inbounds r = (A.chunks[i1] & u) != 0
    return r

end
Base.@propagate_inbounds Base.@assume_effects :noub_if_noinbounds function Base.setindex!(A::FixedSizeBitArray{N}, x, inds::Vararg{Int, N}) where {N}
    @boundscheck checkbounds(A, inds...)
    i1, i2 = get_chunks_id(inds)
    u = UInt64(1) << i2
    @inbounds begin
        c = A.chunks[i1]
        A.chunks[i1] = ifelse(x, c | u, c & ~u)
    end
    return A
end

Base.size(a::FixedSizeBitArray) = a.size
Base.isassigned(a::FixedSizeBitArray, i::Int) = 1 <= i <= length(a)
function Base.fill!(B::FixedSizeBitArray, x)
    y = convert(Bool, x)::Bool
    fill!(B.chunks, UInt64(0)-y)
    # TODO zero partially filled chunks
    return B
end

function (==)(A::FixedSizeBitArray, B::FixedSizeBitArray)
    size(A) != size(B) && return false
    return A.chunks == B.chunks
end
