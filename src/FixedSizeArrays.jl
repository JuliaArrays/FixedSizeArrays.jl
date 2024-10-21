module FixedSizeArrays

export FixedSizeArray, FixedSizeVector, FixedSizeMatrix

public collect_as

"""
    Internal()

Implementation detail. Do not use.
"""
struct Internal end

struct FixedSizeArray{T,N,Mem<:GenericMemory{<:Any,T}} <: DenseArray{T,N}
    mem::Mem
    size::NTuple{N,Int}
    function FixedSizeArray{T,N,M}(::Internal, mem::M, size::NTuple{N,Int}) where {T,N,M<:GenericMemory{<:Any,T}}
        new{T,N,M}(mem, size)
    end
end

function Base.propertynames(
    # the `unused` is here because of https://github.com/JuliaLang/julia/issues/44428
    (@nospecialize unused::FixedSizeArray),
    ::Bool = false,
)
    ()
end

const FixedSizeVector{T} = FixedSizeArray{T,1}
const FixedSizeMatrix{T} = FixedSizeArray{T,2}

const default_underlying_storage_type = Memory

function FixedSizeArray{T,N,V}(::UndefInitializer, size::NTuple{N,Int}) where {T,N,V}
    FixedSizeArray{T,N,V}(Internal(), V(undef, checked_dims(size))::V, size)
end
function FixedSizeArray{T,N,V}(::UndefInitializer, size::NTuple{N,Integer}) where {T,N,V}
    ints = map(Int, size)::NTuple{N,Int}  # prevent infinite recursion
    FixedSizeArray{T,N,V}(undef, ints)
end
function FixedSizeArray{T,N,V}(::UndefInitializer, size::Vararg{Integer,N}) where {T,N,V}
    FixedSizeArray{T,N,V}(undef, size)
end
function FixedSizeArray{T,<:Any,V}(::UndefInitializer, size::NTuple{N,Integer}) where {T,N,V}
    FixedSizeArray{T,N,V}(undef, size)
end
function FixedSizeArray{T,<:Any,V}(::UndefInitializer, size::Vararg{Integer,N}) where {T,N,V}
    FixedSizeArray{T,N,V}(undef, size)
end
function FixedSizeArray{T,N}(::UndefInitializer, size::NTuple{N,Integer}) where {T,N}
    FixedSizeArray{T,N,default_underlying_storage_type{T}}(undef, size)
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
Base.@propagate_inbounds Base.@assume_effects :noub_if_noinbounds function Base.setindex!(A::FixedSizeArray{T}, x, i::Int) where {T}
    @boundscheck checkbounds(A, i)
    @inbounds A.mem[i] = x
    return A
end

Base.size(a::FixedSizeArray) = a.size

function Base.similar(::T, ::Type{E}, size::NTuple{N,Int}) where {T<:FixedSizeArray,E,N}
    with_replaced_parameters(DenseArray, T, Val(E), Val(N))(undef, size)
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

function Base.BroadcastStyle(::Type{T}) where {T<:FixedSizeArray}
    Broadcast.ArrayStyle{stripped_type(DenseArray, T)}()
end

function Base.similar(
    bc::Broadcast.Broadcasted{Broadcast.ArrayStyle{S}},
    ::Type{E},
) where {S<:FixedSizeArray,E}
    similar(S{E}, axes(bc))
end

# helper functions

normalized_type(::Type{T}) where {T} = T

function stripped_type_unchecked(::Type{DenseVector}, ::Type{<:GenericMemory{K,<:Any,AS}}) where {K,AS}
    GenericMemory{K,<:Any,AS}
end

Base.@assume_effects :consistent function stripped_type_unchecked(
    ::Type{DenseArray}, ::Type{<:FixedSizeArray{<:Any,<:Any,V}},
) where {V}
    U = stripped_type(DenseVector, V)
    FixedSizeArray{E,N,U{E}} where {E,N}
end

function stripped_type(::Type{T}, ::Type{S}) where {T,S<:T}
    ret = stripped_type_unchecked(T, S)::Type{<:T}::UnionAll
    S::Type{<:ret}
    normalized_type(ret)  # ensure `UnionAll` type variable order is normalized
end

function with_replaced_parameters(::Type{T}, ::Type{S}, ::Val{P1}, ::Val{P2}) where {T,S<:T,P1,P2}
    t = T{P1,P2}::Type{<:T}
    s = stripped_type(T, S)
    S::Type{<:s}
    s{P1,P2}::Type{<:s}::Type{<:T}::Type{<:t}
end

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
end
function check_count_value(n)
    throw(ArgumentError("count must be an `Int`"))
end

# TODO: use `SpecFSA` for implementing each `FixedSizeArray` constructor?
struct SpecFSA{N,Mem<:GenericMemory} end
function fsa_spec_from_type(::Type{FixedSizeArray})
    SpecFSA{nothing,default_underlying_storage_type}()
end
function fsa_spec_from_type(::Type{FixedSizeArray{<:Any,M}}) where {M}
    check_count_value(M)
    SpecFSA{M,default_underlying_storage_type}()
end
function fsa_spec_from_type(::Type{FixedSizeArray{E}}) where {E}
    E::Type
    SpecFSA{nothing,default_underlying_storage_type{E}}()
end
function fsa_spec_from_type(::Type{FixedSizeArray{E,M}}) where {E,M}
    check_count_value(M)
    E::Type
    SpecFSA{M,default_underlying_storage_type{E}}()
end
function fsa_spec_from_type(::Type{FixedSizeArray{E,<:Any,V}}) where {E,V}
    E::Type
    V::Type{<:DenseVector{E}}
    SpecFSA{nothing,V}()
end
function fsa_spec_from_type(::Type{FixedSizeArray{E,M,V}}) where {E,M,V}
    check_count_value(M)
    E::Type
    V::Type{<:DenseVector{E}}
    SpecFSA{M,V}()
end
for V ∈ (Memory, AtomicMemory)
    T = FixedSizeArray{E,M,V{E}} where {E,M}
    @eval begin
        function fsa_spec_from_type(::Type{$T})
            SpecFSA{nothing,$V}()
        end
        function fsa_spec_from_type(::Type{($T){<:Any,M}}) where {M}
            check_count_value(M)
            SpecFSA{M,$V}()
        end
    end
end

parent_type(::Type{<:FixedSizeArray{<:Any,<:Any,T}}) where {T} = T

underlying_storage(m) = m
underlying_storage(f::FixedSizeArray) = f.mem

axes_are_one_based(axes) = all(isone ∘ first, axes)

# converting constructors for copying other array types

function FixedSizeArray{T,N,V}(src::AbstractArray{S,N}) where {T,N,V,S}
    axs = axes(src)
    if !axes_are_one_based(axs)
        throw(DimensionMismatch("source array has a non-one-based indexing axis"))
    end
    # Can't use `Base.size` because, according to it's doc string, it's not
    # available for all `AbstractArray` types.
    size = map(length, axs)
    dst = FixedSizeArray{T,N,V}(undef, size)
    copyto!(dst.mem, src)
    dst
end

FixedSizeArray{T,<:Any,V}(a::AbstractArray{<:Any,N}) where {V,T,N} = FixedSizeArray{T,N,V}(a)

FixedSizeArray{T,N}(a::AbstractArray{<:Any,N}) where {T,N} = FixedSizeArray{T,N,default_underlying_storage_type{T}}(a)
FixedSizeArray{T}(a::AbstractArray{<:Any,N})   where {T,N} = FixedSizeArray{T,N}(a)
FixedSizeArray{<:Any,N}(a::AbstractArray{T,N}) where {T,N} = FixedSizeArray{T,N}(a)
FixedSizeArray(a::AbstractArray{T,N})          where {T,N} = FixedSizeArray{T,N}(a)

# conversion

Base.convert(::Type{T}, a::T) where {T<:FixedSizeArray} = a
Base.convert(::Type{T}, a::AbstractArray) where {T<:FixedSizeArray} = T(a)::T

# `copyto!`

Base.@propagate_inbounds function copyto5!(dst, doff, src, soff, n)
    copyto!(underlying_storage(dst), doff, underlying_storage(src), soff, n)
    dst
end

Base.@propagate_inbounds function copyto2!(dst, src)
    copyto!(underlying_storage(dst), underlying_storage(src))
    dst
end

Base.@propagate_inbounds Base.copyto!(dst::FixedSizeArray, doff::Integer, src::FixedSizeArray, soff::Integer, n::Integer) = copyto5!(dst, doff, src, soff, n)
Base.@propagate_inbounds Base.copyto!(dst::FixedSizeArray, src::FixedSizeArray) = copyto2!(dst, src)

for A ∈ (Array, GenericMemory)  # Add more? Too bad we have to hardcode to avoid ambiguity.
    @eval begin
        Base.@propagate_inbounds Base.copyto!(dst::FixedSizeArray, doff::Integer, src::$A,             soff::Integer, n::Integer) = copyto5!(dst, doff, src, soff, n)
        Base.@propagate_inbounds Base.copyto!(dst::$A,             doff::Integer, src::FixedSizeArray, soff::Integer, n::Integer) = copyto5!(dst, doff, src, soff, n)

        Base.@propagate_inbounds Base.copyto!(dst::FixedSizeArray, src::$A            ) = copyto2!(dst, src)
        Base.@propagate_inbounds Base.copyto!(dst::$A,             src::FixedSizeArray) = copyto2!(dst, src)
    end
end

# unsafe: the native address of the array's storage

Base.cconvert(::Type{<:Ptr}, a::FixedSizeArray) = a.mem

# `elsize`: part of the strided arrays interface, used for `pointer`

Base.elsize(::Type{A}) where {A<:FixedSizeArray} = Base.elsize(parent_type(A))

# `reshape`: specializing it to ensure it returns a `FixedSizeArray`

function Base.reshape(a::FixedSizeArray{T,<:Any,V}, size::NTuple{N,Int}) where {V,T,N}
    len = checked_dims(size)
    if length(a) != len
        throw(DimensionMismatch("new shape not consistent with existing array length"))
    end
    FixedSizeArray{T,N,V}(Internal(), a.mem, size)
end

# `collect_as`

function collect_as_fsa0(iterator, ::SpecFSA{0,V}) where {V}
    V::UnionAll
    x = only(iterator)
    E = typeof(x)::Type
    ret = FixedSizeArray{E,0,V{E}}(undef)
    ret[] = x
    ret
end

function collect_as_fsa0(iterator, ::SpecFSA{0,V}) where {E,V<:DenseVector{E}}
    E::Type
    x = only(iterator)
    ret = FixedSizeArray{E,0,V}(undef)
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
end

function collect_as_fsam_with_shape(
    iterator, ::SpecFSA{M,V}, shape::Tuple{Vararg{Int}},
) where {M,V}
    V::UnionAll
    E = eltype(iterator)::Type
    U = V{E}
    ret = FixedSizeArray{E,M,U}(undef, shape)
    fill_fsa_from_iterator!(ret, iterator)
    map(identity, ret)::(FixedSizeArray{T,M,V{T}} where {T})
end

function collect_as_fsam_with_shape(
    iterator, ::SpecFSA{M,V}, shape::Tuple{Vararg{Int}},
) where {M,E,V<:DenseVector{E}}
    E::Type
    ret = FixedSizeArray{E,M,V}(undef, shape)
    fill_fsa_from_iterator!(ret, iterator)
    ret::FixedSizeArray{E,M,V}
end

function collect_as_fsam(iterator, spec::SpecFSA{M}) where {M}
    check_count_value(M)
    shape = if isone(M)
        (length(iterator),)
    else
        size(iterator)
    end::NTuple{M,Any}
    shap = map(Int, shape)::NTuple{M,Int}
    collect_as_fsam_with_shape(iterator, spec, shap)::FixedSizeArray{<:Any,M}
end

function collect_as_fsa1_from_unknown_length(iterator, ::SpecFSA{1,V}) where {V}
    V::UnionAll
    v = collect(iterator)::AbstractVector
    t = FixedSizeVector(v)::FixedSizeVector
    s = map(identity, t)::FixedSizeVector  # fix element type
    et = eltype(s)
    FixedSizeVector{et,V{et}}(s)  # fix underlying storage type
end

function collect_as_fsa1_from_unknown_length(iterator, ::SpecFSA{1,V}) where {E,V<:DenseVector{E}}
    E::Type
    v = collect(E, iterator)::AbstractVector{E}
    T = FixedSizeVector{E,V}
    T(v)::T
end

function collect_as_fsa_impl(iterator, spec::SpecFSA{0}, ::LengthIsKnown)
    collect_as_fsa0(iterator, spec)::FixedSizeArray{<:Any,0}
end

function collect_as_fsa_impl(iterator, spec::SpecFSA, ::LengthIsKnown)
    collect_as_fsam(iterator, spec)::FixedSizeArray
end

function collect_as_fsa_impl(iterator, spec::SpecFSA{1}, ::LengthIsUnknown)
    collect_as_fsa1_from_unknown_length(iterator, spec)::FixedSizeVector
end

function collect_as_fsa_checked(iterator, ::SpecFSA{nothing,V}, ::Val{M}, length_status) where {V,M}
    check_count_value(M)
    collect_as_fsa_impl(iterator, SpecFSA{M,V}(), length_status)::FixedSizeArray{<:Any,M}
end

function collect_as_fsa_checked(iterator, spec::SpecFSA{M}, ::Val{M}, length_status) where {M}
    check_count_value(M)
    collect_as_fsa_impl(iterator, spec, length_status)::FixedSizeArray{<:Any,M}
end

"""
    collect_as(t::Type{<:FixedSizeArray}, iterator)

Tries to construct a value of type `t` from the iterator `iterator`. The type `t`
must either be concrete, or a `UnionAll` without constraints.
"""
function collect_as(::Type{T}, iterator) where {T<:FixedSizeArray}
    spec = fsa_spec_from_type(T)::SpecFSA
    size_class = Base.IteratorSize(iterator)
    if size_class == Base.IsInfinite()
        throw(ArgumentError("iterator is infinite, can't fit infinitely many elements into a `FixedSizeArray`"))
    end
    dim_count_int = dimension_count_of(size_class)
    check_count_value(dim_count_int)
    dim_count = Val(dim_count_int)::Val
    len_stat = length_status(size_class)
    collect_as_fsa_checked(iterator, spec, dim_count, len_stat)::T
end

if isdefined(Base, :dataids) && (Base.dataids isa Function)
    # This is an internal, non-public function which is nevertheless needed to
    # get good performance in some cases (e.g. broadcasting):
    # <https://github.com/JuliaArrays/FixedSizeArrays.jl/issues/63>.
    Base.dataids(a::FixedSizeArray) = Base.dataids(a.mem)
end

end # module FixedSizeArrays
