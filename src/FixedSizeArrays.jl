module FixedSizeArrays

export FixedSizeArray, FixedSizeVector, FixedSizeMatrix

struct FixedSizeArray{T,N} <: DenseArray{T,N}
    mem::Memory{T}
    size::NTuple{N,Int}
    function FixedSizeArray{T,N}(mem::Memory{T}, size::NTuple{N,Int}) where {T,N}
        new{T,N}(mem, size)
    end
end

const FixedSizeVector{T} = FixedSizeArray{T,1}
const FixedSizeMatrix{T} = FixedSizeArray{T,2}

function FixedSizeArray{T,N}(::UndefInitializer, size::NTuple{N,Int}) where {T,N}
    FixedSizeArray{T,N}(Memory{T}(undef, Core.checked_dims(size...)), size)
end
function FixedSizeArray{T,N}(::UndefInitializer, size::Vararg{Int,N}) where {T,N}
    FixedSizeArray{T,N}(undef, size)
end
function FixedSizeArray{T}(::UndefInitializer, size::Vararg{Int,N}) where {T,N}
    FixedSizeArray{T,N}(undef, size)
end
function FixedSizeArray{T}(::UndefInitializer, size::NTuple{N,Int}) where {T,N}
    FixedSizeArray{T,N}(undef, size)
end

Base.IndexStyle(::Type{<:FixedSizeArray}) = IndexLinear()
Base.@propagate_inbounds Base.getindex(A::FixedSizeArray, i::Int) = A.mem[i]
Base.@propagate_inbounds Base.setindex!(A::FixedSizeArray, v, i::Int) = A.mem[i] = v

Base.size(a::FixedSizeArray) = getfield(a, :size)

function Base.similar(::FixedSizeArray, ::Type{S}, size::NTuple{N,Int}) where {S,N}
    FixedSizeArray{S,N}(undef, size...)
end

Base.isassigned(a::FixedSizeArray, i::Int) = isassigned(a.mem, i)

# broadcasting

function Base.BroadcastStyle(::Type{<:FixedSizeArray})
    Broadcast.ArrayStyle{FixedSizeArray}()
end

function Base.similar(
    bc::Broadcast.Broadcasted{Broadcast.ArrayStyle{FixedSizeArray}},
    ::Type{E},
) where {E}
    similar(FixedSizeArray{E}, axes(bc))
end

# helper functions

memory_of(m::Memory) = m
memory_of(f::FixedSizeArray) = f.mem

first_linear_index(a) = first(eachindex(IndexLinear(), a))

axes_are_one_based(axes) = all(isone ∘ first, axes)

# converting constructors for copying other array types

function FixedSizeArray{T,N}(src::AbstractArray{S,N}) where {T,N,S}
    axs = axes(src)
    axes_are_one_based(axs) ||
        throw(DimensionMismatch("source array has a non-one-based indexing axis"))
    # Can't use `Base.size` because, according to it's doc string, it's not
    # available for all `AbstractArray` types.
    size = map(length, axs)
    dst = FixedSizeArray{T,N}(undef, size)
    copyto!(dst, src)::FixedSizeArray{T,N}
end

FixedSizeArray{T}(a::AbstractArray{<:Any,N})   where {T,N} = FixedSizeArray{T,N}(a)
FixedSizeArray{<:Any,N}(a::AbstractArray{T,N}) where {T,N} = FixedSizeArray{T,N}(a)
FixedSizeArray(a::AbstractArray{T,N})          where {T,N} = FixedSizeArray{T,N}(a)

# conversion

Base.convert(::Type{T}, a::T) where {T<:FixedSizeArray} = a
Base.convert(::Type{T}, a::AbstractArray) where {T<:FixedSizeArray} = T(a)::T

# `copyto!` between `FixedSizeArray` and `Memory`

Base.@propagate_inbounds function copyto5!(dst, doff, src, soff, n)
    if !iszero(n)
        (n < false) &&
            throw(ArgumentError("the number of elements to copy must be nonnegative"))
        @boundscheck checkbounds(dst, doff:doff+n-1)
        @boundscheck checkbounds(src, soff:soff+n-1)
        @inbounds let d, s
            d = GenericMemoryRef(memory_of(dst), doff)
            s = GenericMemoryRef(memory_of(src), soff)
            unsafe_copyto!(d, s, n)
        end
    end
    dst
end

Base.@propagate_inbounds function copyto2!(dst::T, src) where {T}
    fd = first_linear_index(dst)
    fs = first_linear_index(src)
    len = length(src)
    copyto5!(dst, fd, src, fs, len)::T
end

Base.@propagate_inbounds Base.copyto!(dst::FixedSizeArray, doff::Integer, src::FixedSizeArray, soff::Integer, n::Integer) = copyto5!(dst, doff, src, soff, n)
Base.@propagate_inbounds Base.copyto!(dst::FixedSizeArray, doff::Integer, src::Memory        , soff::Integer, n::Integer) = copyto5!(dst, doff, src, soff, n)
Base.@propagate_inbounds Base.copyto!(dst::Memory        , doff::Integer, src::FixedSizeArray, soff::Integer, n::Integer) = copyto5!(dst, doff, src, soff, n)

Base.@propagate_inbounds Base.copyto!(dst::FixedSizeArray, src::FixedSizeArray) = copyto2!(dst, src)
Base.@propagate_inbounds Base.copyto!(dst::FixedSizeArray, src::Memory        ) = copyto2!(dst, src)
Base.@propagate_inbounds Base.copyto!(dst::Memory        , src::FixedSizeArray) = copyto2!(dst, src)

# unsafe: the native address of the array's storage

Base.unsafe_convert(::Type{Ptr{T}}, a::FixedSizeArray{T}) where {T} = Base.unsafe_convert(Ptr{T}, a.mem)

end # module FixedSizeArrays
