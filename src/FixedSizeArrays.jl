module FixedSizeArrays

export FixedSizeArray, FixedSizeVector, FixedSizeMatrix

"""
    Internal()

Implementation detail. Do not use.
"""
struct Internal end

struct FixedSizeArray{T,N} <: DenseArray{T,N}
    mem::Memory{T}
    size::NTuple{N,Int}
    function FixedSizeArray{T,N}(::Internal, mem::Memory{T}, size::NTuple{N,Int}) where {T,N}
        new{T,N}(mem, size)
    end
end

const FixedSizeVector{T} = FixedSizeArray{T,1}
const FixedSizeMatrix{T} = FixedSizeArray{T,2}

function FixedSizeArray{T,N}(::UndefInitializer, size::NTuple{N,Int}) where {T,N}
    FixedSizeArray{T,N}(Internal(), Memory{T}(undef, checked_dims(size)), size)
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

# safe product of a tuple of integers, for calculating dimensions size

checked_dims_impl(a::Int, ::Tuple{}, have_overflow::Bool) = (a, have_overflow)
function checked_dims_impl(a::Int, t::Tuple{Int,Vararg{Int,N}}, have_overflow::Bool) where {N}
    b = first(t)
    (m, o) = Base.Checked.mul_with_overflow(a, b)
    r = Base.tail(t)::NTuple{N,Int}
    checked_dims_impl(m, r, have_overflow | o)::Tuple{Int,Bool}
end

checked_dims(::Tuple{}) = 1
function checked_dims(t::Tuple{Int,Vararg{Int,N}}) where {N}
    any_is_zero     = any(iszero,                 t)::Bool
    any_is_negative = any((x -> x < false),       t)::Bool
    any_is_typemax  = any((x -> x == typemax(x)), t)::Bool
    a = first(t)
    r = Base.tail(t)::NTuple{N,Int}
    (product, have_overflow) = checked_dims_impl(a, r, false)::Tuple{Int,Bool}
    if any_is_negative
        throw(ArgumentError("array dimension size can't be negative"))
    end
    if any_is_typemax
        throw(ArgumentError("array dimension size can't be the maximum representable value"))
    end
    if have_overflow & !any_is_zero
        throw(ArgumentError("array dimensions too great, can't represent length"))
    end
    product
end

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

parent_type(::Type{<:FixedSizeArray{T}}) where {T} = Memory{T}

memory_of(m::Memory) = m
memory_of(f::FixedSizeArray) = f.mem

first_linear_index(a) = first(eachindex(IndexLinear(), a))

axes_are_one_based(axes) = all(isone ∘ first, axes)

# converting constructors for copying other array types

function FixedSizeArray{T,N}(src::AbstractArray{S,N}) where {T,N,S}
    axs = axes(src)
    if !axes_are_one_based(axs)
        throw(DimensionMismatch("source array has a non-one-based indexing axis"))
    end
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
        if n < false
            throw(ArgumentError("the number of elements to copy must be nonnegative"))
        end
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

# `elsize`: part of the strided arrays interface, used for `pointer`

Base.elsize(::Type{A}) where {A<:FixedSizeArray} = Base.elsize(parent_type(A))

# `reshape`: specializing it to ensure it returns a `FixedSizeArray`

function Base.reshape(a::FixedSizeArray{T}, size::NTuple{N,Int}) where {T,N}
    len = checked_dims(size)
    if length(a) != len
        throw(DimensionMismatch("new shape not consistent with existing array length"))
    end
    FixedSizeArray{T,N}(Internal(), a.mem, size)
end

end # module FixedSizeArrays
