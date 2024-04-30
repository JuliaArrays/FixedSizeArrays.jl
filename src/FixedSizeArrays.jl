module FixedSizeArrays

export FixedSizeArray, FixedSizeVector, FixedSizeMatrix

public collect_as

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
function FixedSizeArray{T,N}(::UndefInitializer, size::NTuple{N,Integer}) where {T,N}
    ints = map(Int, size)::NTuple{N,Int}  # prevent infinite recursion
    FixedSizeArray{T,N}(undef, ints)
end
function FixedSizeArray{T,N}(::UndefInitializer, size::Vararg{Integer,N}) where {T,N}
    FixedSizeArray{T,N}(undef, size)
end
function FixedSizeArray{T}(::UndefInitializer, size::Vararg{Integer,N}) where {T,N}
    FixedSizeArray{T,N}(undef, size)
end
function FixedSizeArray{T}(::UndefInitializer, size::NTuple{N,Integer}) where {T,N}
    FixedSizeArray{T,N}(undef, size)
end

Base.IndexStyle(::Type{<:FixedSizeArray}) = IndexLinear()
Base.@propagate_inbounds Base.getindex(A::FixedSizeArray, i::Int) = A.mem[i]
Base.@propagate_inbounds Base.setindex!(A::FixedSizeArray, v, i::Int) = A.mem[i] = v

Base.size(a::FixedSizeArray) = a.size

function Base.similar(::FixedSizeArray, ::Type{S}, size::NTuple{N,Int}) where {S,N}
    FixedSizeArray{S,N}(undef, size)
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

dimension_count_of(::Base.SizeUnknown) = 1
dimension_count_of(::Base.HasLength) = 1
dimension_count_of(::Base.HasShape{N}) where {N} = convert(Int, N)::Int

struct LengthIsUnknown end
struct LengthIsKnown end
length_status(::Base.SizeUnknown) = LengthIsUnknown()
length_status(::Base.HasLength) = LengthIsKnown()
length_status(::Base.HasShape) = LengthIsKnown()

function check_count_value(n::Int)
    if n < 0
        throw(ArgumentError("count can't be negative"))
    end
    nothing
end
function check_count_value(n)
    throw(ArgumentError("count must be an `Int`"))
end

struct SpecFSA{T,N} end
function fsa_spec_to_type(::SpecFSA{nothing,nothing})
    FixedSizeArray
end
function fsa_spec_to_type(::SpecFSA{nothing,M}) where {M}
    check_count_value(M)
    FixedSizeArray{<:Any,M}
end
function fsa_spec_to_type(::SpecFSA{E,nothing}) where {E}
    FixedSizeArray{E::Type}
end
function fsa_spec_to_type(::SpecFSA{E,M}) where {E,M}
    check_count_value(M)
    FixedSizeArray{E::Type,M}
end
function fsa_spec_from_type(::Type{FixedSizeArray})
    SpecFSA{nothing,nothing}()
end
function fsa_spec_from_type(::Type{FixedSizeArray{<:Any,M}}) where {M}
    check_count_value(M)
    SpecFSA{nothing,M}()
end
function fsa_spec_from_type(::Type{FixedSizeArray{E}}) where {E}
    SpecFSA{E::Type,nothing}()
end
function fsa_spec_from_type(::Type{FixedSizeArray{E,M}}) where {E,M}
    check_count_value(M)
    SpecFSA{E::Type,M}()
end

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

# `collect_as`

function collect_as_fsa0(iterator, ::Val{nothing})
    x = only(iterator)
    ret = FixedSizeArray{typeof(x),0}(undef)
    ret[] = x
    ret
end

function collect_as_fsa0(iterator, ::Val{E}) where {E}
    E::Type
    x = only(iterator)
    ret = FixedSizeArray{E,0}(undef)
    ret[] = x
    ret
end

function fill_fsa_from_iterator!(a, iterator)
    actual_count = 0
    for e ∈ iterator
        actual_count += 1
        a[actual_count] = e
    end
    if actual_count != length(a)
        throw(ArgumentError("`size`-`length` inconsistency"))
    end
    nothing
end

function collect_as_fsam_with_shape(
    iterator, ::SpecFSA{nothing,M}, shape::Tuple{Vararg{Int}},
) where {M}
    E = eltype(iterator)::Type
    ret = FixedSizeArray{E,M}(undef, shape)
    fill_fsa_from_iterator!(ret, iterator)
    map(identity, ret)::FixedSizeArray{<:Any,M}
end

function collect_as_fsam_with_shape(
    iterator, ::SpecFSA{E,M}, shape::Tuple{Vararg{Int}},
) where {E,M}
    E::Type
    ret = FixedSizeArray{E,M}(undef, shape)
    fill_fsa_from_iterator!(ret, iterator)
    ret::FixedSizeArray{E,M}
end

function collect_as_fsam(iterator, spec::SpecFSA{<:Any,M}) where {M}
    check_count_value(M)
    shape = if isone(M)
        (length(iterator),)
    else
        size(iterator)
    end::NTuple{M,Any}
    shap = map(Int, shape)::NTuple{M,Int}
    collect_as_fsam_with_shape(iterator, spec, shap)::FixedSizeArray{<:Any,M}
end

function collect_as_fsa1_from_unknown_length(iterator, ::Val{nothing})
    v = collect(iterator)::AbstractVector
    T = FixedSizeVector
    map(identity, T(v))::T
end

function collect_as_fsa1_from_unknown_length(iterator, ::Val{E}) where {E}
    E::Type
    v = collect(E, iterator)::AbstractVector{E}
    T = FixedSizeVector{E}
    T(v)::T
end

function collect_as_fsa_impl(iterator, ::SpecFSA{E,0}, ::LengthIsKnown) where {E}
    collect_as_fsa0(iterator, Val(E))::FixedSizeArray{<:Any,0}
end

function collect_as_fsa_impl(iterator, spec::SpecFSA, ::LengthIsKnown)
    collect_as_fsam(iterator, spec)::FixedSizeArray
end

function collect_as_fsa_impl(iterator, ::SpecFSA{E,1}, ::LengthIsUnknown) where {E}
    collect_as_fsa1_from_unknown_length(iterator, Val(E))::FixedSizeVector
end

function collect_as_fsa_checked(iterator, ::SpecFSA{E,nothing}, ::Val{M}, length_status) where {E,M}
    check_count_value(M)
    collect_as_fsa_impl(iterator, SpecFSA{E,M}(), length_status)::FixedSizeArray{<:Any,M}
end

function collect_as_fsa_checked(iterator, ::SpecFSA{E,M}, ::Val{M}, length_status) where {E,M}
    check_count_value(M)
    collect_as_fsa_impl(iterator, SpecFSA{E,M}(), length_status)::FixedSizeArray{<:Any,M}
end

function collect_as_fsa(iterator, spec::SpecFSA)
    size_class = Base.IteratorSize(iterator)
    if size_class == Base.IsInfinite()
        throw(ArgumentError("iterator is infinite, can't fit infinitely many elements into a `FixedSizeArray`"))
    end
    dim_count_int = dimension_count_of(size_class)
    check_count_value(dim_count_int)
    dim_count = Val(dim_count_int)::Val
    len_stat = length_status(size_class)
    R = fsa_spec_to_type(spec)::Type{<:FixedSizeArray}
    collect_as_fsa_checked(iterator, spec, dim_count, len_stat)::R
end

"""
    collect_as(t::Type{<:FixedSizeArray}, iterator)

Tries to construct a value of type `t` from the iterator `iterator`. The type `t`
must either be concrete, or a `UnionAll` without constraints.
"""
function collect_as(::Type{T}, iterator) where {T<:FixedSizeArray}
    spec = fsa_spec_from_type(T)::SpecFSA
    collect_as_fsa(iterator, spec)::T
end

end # module FixedSizeArrays
