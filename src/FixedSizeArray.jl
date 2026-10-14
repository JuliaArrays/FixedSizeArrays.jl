export
    BoundsErrorLight,
    FixedSizeArray, FixedSizeVector, FixedSizeMatrix,
    FixedSizeArrayDefault, FixedSizeVectorDefault, FixedSizeMatrixDefault

"""
    FixedSizeArray{T,N,Mem<:DenseVector{T}}(undef, size::NTuple{N,Int})
    FixedSizeArray{T,N,Mem<:DenseVector{T}}(undef, size1::Int, size2::Int, ...)
    FixedSizeArray{T,N,Mem<:DenseVector{T}}(array::AbstractArray)

Construct a fixed-size `FixedSizeArray` with element type `T`, number of dimensions `N`, and memory backend `Mem`.

The memory backend `Mem` is a `DenseVector{T}` (e.g. `Memory{T}` or `Vector{T}`) of length
`prod(size)`, or, on Julia v1.11 and later, a `GenericMemoryRef{kind,T}` (e.g. `MemoryRef{T}`)
to the first element of the array, with at least `prod(size)` elements from it on in its memory.

Convenient aliases provided:

* [`FixedSizeVector{T,Mem}`](@ref): for 1-dimensional fixed-size arrays
* [`FixedSizeMatrix{T,Mem}`](@ref): for 2-dimensional fixed-size arrays
* [`FixedSizeArrayDefault{T,N}`](@ref): for `N`-dimensional fixed-size arrays, with default memory backend
* [`FixedSizeVectorDefault{T}`](@ref): for 1-dimensional fixed-size arrays, with default memory backend
* [`FixedSizeMatrixDefault{T}`](@ref): for 2-dimensional fixed-size arrays, with default memory backend
"""
struct FixedSizeArray{T,N,Mem<:StorageBound{T}} <: DenseArray{T,N}
    mem::Mem
    size::NTuple{N,Int}
    global function new_fixed_size_array(mem::StorageBound{T}, size::NTuple{N,Int}) where {T,N}
        new{T,N,typeof(mem)}(mem, size)
    end
end

function Base.propertynames(
    # the `unused` is here because of https://github.com/JuliaLang/julia/issues/44428
    (@nospecialize unused::FixedSizeArray),
    ::Bool = false,
)
    ()
end

const FixedSizeArrayAllowedConstructorType = let
    function f(Mem::UnionAll)
        Union{
            # 0 fixed parameters
            Type{FixedSizeArray{T, N, Mem{T}} where {T, N}},
            # 1 fixed parameter
            Type{FixedSizeArray{T, N, Mem{T}} where {T}} where {N},
        }
    end
    special_storage_types = (Vector, optional_memory..., optional_memory_ref...)
    Union{
        # 0 fixed parameters
        Type{FixedSizeArray},
        # 1 fixed parameter
        Type{FixedSizeArray{T}} where {T},
        Type{FixedSizeArray{T, N} where {T}} where {N},
        # 2 fixed parameters
        Type{FixedSizeArray{T, N}} where {T, N},
        Type{FixedSizeArray{T, N, Mem} where {N}} where {T, Mem <: StorageBound{T}},
        # 3 fixed parameters
        Type{FixedSizeArray{T, N, Mem}} where {T, N, Mem <: StorageBound{T}},
        # special cases depending on the underlying storage type
        map(f, special_storage_types)...,
    }
end

function check_constructor_is_allowed(::Type{T}) where {T <: FixedSizeArray}
    if Type{T} <: FixedSizeArrayAllowedConstructorType
        T
    else
        throw(ArgumentError("technical limitation: type is not among the allowed `FixedSizeArray` constructors (try deleting type parameter constraints)"))
    end
end

"""
    FixedSizeVector{T,Mem<:DenseVector{T}}(undef, size::NTuple{1,Int})
    FixedSizeVector{T,Mem<:DenseVector{T}}(undef, size1::Int)
    FixedSizeVector{T,Mem<:DenseVector{T}}(array::AbstractVector)

Construct a fixed-size 1-dimensional [`FixedSizeArray`](@ref) with element type `T`, and memory backend `Mem`.

Convenient aliases provided:

* [`FixedSizeVectorDefault{T}`](@ref): for 1-dimensional fixed-size arrays, with default memory backend
"""
const FixedSizeVector{T} = FixedSizeArray{T,1}
"""
    FixedSizeMatrix{T,Mem<:DenseMatrix{T}}(undef, size::NTuple{2,Int})
    FixedSizeMatrix{T,Mem<:DenseMatrix{T}}(undef, size1::Int, size2::Int)
    FixedSizeMatrix{T,Mem<:DenseMatrix{T}}(array::AbstractMatrix)

Construct a fixed-size 2-dimensional [`FixedSizeArray`](@ref) with element type `T`, and memory backend `Mem`.

Convenient aliases provided:

* [`FixedSizeMatrixDefault{T}`](@ref): for 2-dimensional fixed-size arrays, with default memory backend
"""
const FixedSizeMatrix{T} = FixedSizeArray{T,2}

function parent_type_with_default(::Type{<:(FixedSizeArray{T, N, Mem} where {N})}) where {T, Mem <: StorageBound{T}}
    Mem
end
function parent_type_with_default(::Type{<:FixedSizeArray{T}}) where {T}
    default_underlying_storage_type{T}
end
function parent_type_with_default(::Type{<:FixedSizeArray})
    default_underlying_storage_type
end
for Mem ∈ (Vector, optional_memory..., optional_memory_ref...)
    FSA = FixedSizeArray{T, N, Mem{T}} where {T, N}
    t_fsa = Union{
        Type{FSA},
        Type{FSA{T, N} where {T}} where {N},
    }
    @eval begin
        function parent_type_with_default(::$t_fsa)
            $Mem
        end
    end
end

function check_ndims(::Type{<:(FixedSizeArray{T, N} where {T})}, size::Tuple{Vararg{Integer}}) where {N}
    if size isa Tuple{Vararg{Any, N}}
        size
    else
        throw(DimensionMismatch("mismatch between dimension count in type and the length of the size tuple"))
    end
end
function check_ndims(::Type{<:FixedSizeArray}, size::Tuple{Vararg{Integer}})
    size
end

function undef_constructor(::Type{FSA}, size::Tuple{Vararg{Integer}}) where {T, FSA <: FixedSizeArray{T}}
    size = check_ndims(FSA, size)
    s = map(Int, size)
    Mem = parent_type_with_default(FSA)
    len = checked_dims(s)
    mem = new_storage(Mem, len)
    new_fixed_size_array(mem, s)
end

function (::Type{FSA})(::UndefInitializer, size::Tuple{Vararg{Integer}}) where {T, FSA <: FixedSizeArray{T}}
    undef_constructor(FSA, size)
end
function (::Type{FSA})(::UndefInitializer, size::Vararg{Integer}) where {T, FSA <: FixedSizeArray{T}}
    undef_constructor(FSA, size)
end

macro assume_noub_if_noinbounds(x)
    expr = :(
        let f
            Base.@assume_effects :noub_if_noinbounds f() = 7
        end
    )
    effect_is_recognized = try
        eval(expr)
        true
    catch
        false
    end
    if effect_is_recognized
        :(Base.@assume_effects :noub_if_noinbounds $x)
    else
        x
    end
end

Base.IndexStyle(::Type{<:FixedSizeArray}) = IndexLinear()
Base.@propagate_inbounds function Base.getindex(A::FixedSizeArray, i::Int)
    @boundscheck checkbounds(A, i)
    @inbounds storage_getindex(parent(A), i)
end
Base.@propagate_inbounds @assume_noub_if_noinbounds function Base.setindex!(A::FixedSizeArray, x, i::Int)
    @boundscheck checkbounds(A, i)
    @inbounds storage_setindex!(parent(A), x, i)
    return A
end

# Access to the memory backend. The bounds are checked by the `FixedSizeArray`, these
# methods are always called in an `@inbounds` context. Methods for `GenericMemoryRef` are
# in `memoryref.jl`.
Base.@propagate_inbounds storage_getindex(m, i::Int) = m[i]
Base.@propagate_inbounds function storage_setindex!(m, x, i::Int)
    m[i] = x
    nothing
end
Base.@propagate_inbounds storage_isassigned(m, i::Int) = isassigned(m, i)
new_storage(::Type{Mem}, len::Int) where {Mem} = Mem(undef, len)
copy_storage(m, ::Int) = copy(m)
collect_storage(collect, ::Type{Mem}, iterator) where {Mem} = collect(Mem, iterator)
storage_elsize(::Type{Mem}) where {Mem} = Base.elsize(Mem)
storage_length(m) = length(m)

function Base.checkbounds(A::FixedSizeArray, is...)
    checkbounds_lightboundserror(A, is...)
end

# Load all dimensions at once. With a plain `getfield(a, :size)`, Julia's codegen loads
# each dimension of a by-reference `FixedSizeArray` argument lazily, at its use. In a
# short-circuiting check such as `axes(a) == expected || throw(...)` the loads of
# dimensions 2:N then only happen in the conditional blocks of the `&&` chain, so they do
# not dominate later uses in loops. LLVM's LICM later hoists the loads of the loop out of
# conditionally executed inner loops and drops their `!tbaa`/`!invariant.load` metadata,
# after which GVN can neither CSE them with the check nor use the checked values, so the
# dimensions stay dynamic in the loop. Passing each dimension through an intrinsic forces
# the tuple to be materialized here (as for `Array` and `PtrArray`) and is a no-op otherwise.
Base.size(a::FixedSizeArray) = map(materialize_int, getfield(a, :size))
materialize_int(n::Int) = Core.Intrinsics.bitcast(Int, n)

function Base.similar(::T, ::Type{E}, size::NTuple{N,Int}) where {T<:FixedSizeArray,E,N}
    spec = TypeParametersElementTypeAndDimensionality()
    S = val_parameter(with_stripped_type_parameters(spec, T)){E, N}
    S(undef, size)
end

function Base.isassigned(a::FixedSizeArray, i::Int)
    @boundscheck checkbounds(Bool, a, i) || return false
    @inbounds storage_isassigned(parent(a), i)
end

# safe product of a tuple of integers, for calculating dimensions size

@noinline function throw_checked_dims(shape::Tuple)
    msg = lazy"""
    invalid array dimensions $shape, check if:
    * any dim is negative
    * any dim is `typemax(T)`
    * trying to compute the total length overflows"""
    throw(ArgumentError(msg))
end

checked_dims(::Tuple{}) = 1
function checked_dims(t::Tuple{Int,Vararg{Int,N}}) where {N}
    product = checked_size_product(t)
    if product === nothing
        throw_checked_dims(t)
    end
    product
end

# broadcasting

struct FixedSizeArrayBroadcastStyle{N, Mem <: StorageBoundAny} <: Broadcast.AbstractArrayStyle{N} end

function (::Type{<:(FixedSizeArrayBroadcastStyle{N, Mem} where {N})})(::Val{M}) where {Mem, M}
    m = check_count_value(M)
    FixedSizeArrayBroadcastStyle{m, Mem}()  # `FixedSizeArray` supports any dimensionality
end

function Base.BroadcastStyle(::Union{Type{<:(FixedSizeArray{T, N, Mem} where {T})},
                                     Type{<:(SubArray{T, N, FixedSizeArray{T, N, Mem}} where {T})}}) where {N, Mem}
    n = check_count_value(N)
    spec = TypeParametersElementType()
    mem = val_parameter(with_stripped_type_parameters(spec, Mem))
    FixedSizeArrayBroadcastStyle{n, mem}()
end

function Base.similar(
    bc::Broadcast.Broadcasted{FixedSizeArrayBroadcastStyle{N, Mem}},
    ::Type{T},
) where {N,Mem,T}
    similar(FixedSizeArray{T, N, Mem{T}}, axes(bc))
end

# helper functions

val_parameter(::Val{P}) where {P} = P

struct TypeParametersElementType end
struct TypeParametersElementTypeAndDimensionality end

"""
    with_stripped_type_parameters_unchecked(spec, t::Type)::Val{s}

An implementation detail of [`with_stripped_type_parameters`](@ref). Don't call
directly.
"""
function with_stripped_type_parameters_unchecked end
function with_stripped_type_parameters_unchecked(::TypeParametersElementType, ::Type{<:Vector})
    s = Vector
    Val{s}()
end
(@isdefined GenericMemory) &&
function with_stripped_type_parameters_unchecked(::TypeParametersElementType, ::Type{<:(GenericMemory{K, T, AS} where {T})}) where {K, AS}
    s = GenericMemory{K, T, AS} where {T}
    Val{s}()
end

# `Base.@assume_effects :consistent` is a workaround for:
# https://github.com/JuliaLang/julia/issues/56966
Base.@assume_effects :consistent function with_stripped_type_parameters_unchecked(::TypeParametersElementTypeAndDimensionality, ::Type{<:(FixedSizeArray{T, N, Mem} where {T, N})}) where {Mem}
    spec_mem = TypeParametersElementType()
    mem_v = with_stripped_type_parameters(spec_mem, Mem)
    mem = val_parameter(mem_v)
    s = FixedSizeArray{T, N, mem{T}} where {T, N}
    Val{s}()
end

"""
    with_stripped_type_parameters(spec, t::Type)::Val{s}

The type `s` is a `UnionAll` supertype of `t`:

```julia
(s isa UnionAll) && (t <: s)
```

Furthermore, `s` has type variables in place of the type parameters specified
via `spec`.

NB: `Val{s}()` is returned instead of `s` so the method would be *consistent*
from the point of view of Julia's effect inference, enabling constant folding.

NB: this function is supposed to only have the one method. To add
functionality, add methods to [`with_stripped_type_parameters_unchecked`](@ref).
"""
function with_stripped_type_parameters(spec, t::Type)
    ret = with_stripped_type_parameters_unchecked(spec, t)
    s = val_parameter(ret)
    s = s::UnionAll
    Val{s}()
end

dimension_count_of(::Base.SizeUnknown) = 1
dimension_count_of(::Base.HasLength) = 1
dimension_count_of(::Base.HasShape{N}) where {N} = convert(Int, N)::Int

function check_count_value(n)
    n = n::Int
    if n < 0
        throw(ArgumentError("count can't be negative"))
    end
    n
end

parent_type(::Type{<:FixedSizeArray{<:Any,<:Any,Mem}}) where {Mem} = Mem

"""
    parent(f::FixedSizeArray)

Return the underlying parent object of the fixed-size array `f`, i.e. its memory backend.
Use `parentindices(f)` to get the valid indices of the parent of `f`.

For a `GenericMemoryRef` backend this is the memory reference to the first element of `f`
(which is not an `AbstractArray`), and `parentindices(f)` are the indices
`(Base.OneTo(length(f)),)` relative to it.
"""
Base.parent(f::FixedSizeArray) = getfield(f, :mem)

Base.parentindices(f::FixedSizeArray) = storage_parentindices(parent(f), f)
storage_parentindices(m, ::FixedSizeArray) = axes(m)

"""
    underlying_storage(f)

The memory backend of `f` as an `AbstractVector` of length `length(f)`, sharing the
elements of `f`: `parent(f)` for a `DenseVector` backend, a view of the memory for a
`GenericMemoryRef` backend.
"""
underlying_storage(m) = m
underlying_storage(f::FixedSizeArray) = storage_vector(parent(f), length(f))
storage_vector(m, ::Int) = m

axes_are_one_based(axes) = all(isone ∘ first, axes)

# converting constructors for copying other array types

function converting_constructor(::Type{FSA}, src::AbstractArray) where {FSA <: FixedSizeArray}
    axs = axes(src)
    if !axes_are_one_based(axs)
        throw(DimensionMismatch("source array has a non-one-based indexing axis"))
    end
    collect_as_haseltype(FSA, src)
end

function (::Type{FSA})(src::AbstractArray) where {N, FSA <: FixedSizeArray{<:Any, N}}
    if N::Int != ndims(src)
        throw(DimensionMismatch("dimensionality mismatch"))
    end
    converting_constructor(FSA, src)
end

function (::Type{FSA})(src::AbstractArray) where {FSA <: FixedSizeArray}
    converting_constructor(FSA, src)
end

# `copy`: avoid the `similar` and `copyto!` in the generic fallback

function Base.copy(a::FixedSizeArray)
    new_fixed_size_array(copy_storage(Base.parent(a), length(a)), size(a))
end

# `getindex` with a `Colon` for the index: better effects than with the generic fallback

function Base.getindex(a::FixedSizeArray, ::Colon)
    c = copy(a)
    if a isa AbstractVector
        c
    else
        reshape(c, :)
    end
end

# `copyto!`

Base.@propagate_inbounds function copyto5!(dst, doff, src, soff, n)
    if use_memoryref_copy(dst, src)
        copyto_memoryref!(dst, Int(doff), src, Int(soff), Int(n))
    else
        copyto!(underlying_storage(dst), doff, underlying_storage(src), soff, n)
    end
    dst
end

Base.@propagate_inbounds function copyto2!(dst, src)
    if use_memoryref_copy(dst, src)
        copyto_memoryref!(dst, 1, src, 1, length(src))
    else
        copyto!(underlying_storage(dst), underlying_storage(src))
    end
    dst
end

# whether to copy directly between memory references, see `memoryref.jl`
use_memoryref_copy(dst, src) = false

Base.@propagate_inbounds Base.copyto!(dst::FixedSizeArray, doff::Integer, src::FixedSizeArray, soff::Integer, n::Integer) = copyto5!(dst, doff, src, soff, n)
Base.@propagate_inbounds Base.copyto!(dst::FixedSizeArray, src::FixedSizeArray) = copyto2!(dst, src)

for A ∈ (Array, optional_generic_memory...)  # Add more? Too bad we have to hardcode to avoid ambiguity.
    @eval begin
        Base.@propagate_inbounds Base.copyto!(dst::FixedSizeArray, doff::Integer, src::$A,             soff::Integer, n::Integer) = copyto5!(dst, doff, src, soff, n)
        Base.@propagate_inbounds Base.copyto!(dst::$A,             doff::Integer, src::FixedSizeArray, soff::Integer, n::Integer) = copyto5!(dst, doff, src, soff, n)

        Base.@propagate_inbounds Base.copyto!(dst::FixedSizeArray, src::$A            ) = copyto2!(dst, src)
        Base.@propagate_inbounds Base.copyto!(dst::$A,             src::FixedSizeArray) = copyto2!(dst, src)
    end
end

# `copy!`: avoid the `resize!` branch of `copy!(::AbstractVector, ::AbstractVector)`

function Base.copy!(dst::FixedSizeVector, src::AbstractVector)
    axes(dst) == axes(src) || throw(ArgumentError(lazy"arrays must have the same axes for copy! (consider using `copyto!`), got $(axes(dst)) and $(axes(src))"))
    copyto!(dst, src)
end

# unsafe: the native address of the array's storage

Base.cconvert(::Type{<:Ptr}, a::FixedSizeArray) = parent(a)

function Base.unsafe_convert(::Type{Ptr{T}}, a::FixedSizeArray{T}) where {T}
    Base.unsafe_convert(Ptr{T}, parent(a))
end

# `elsize`: part of the strided arrays interface, used for `pointer`

Base.elsize(::Type{A}) where {A<:FixedSizeArray} = storage_elsize(parent_type(A))

# `reshape`: specializing it to ensure it returns a `FixedSizeArray`

function Base.reshape(a::FixedSizeArray, size::(NTuple{N,Int} where {N}))
    len = checked_dims(size)
    if length(a) != len
        throw(DimensionMismatch("new shape not consistent with existing array length"))
    end
    new_fixed_size_array(parent(a), size)
end

# `wrap`: construct a `FixedSizeArray` around an existing memory backend, without copying

"""
    FixedSizeArrays.wrap(FixedSizeArray, parent::DenseVector, size::NTuple{N,Integer})
    FixedSizeArrays.wrap(FixedSizeArray, parent::DenseVector, size1::Integer, size2::Integer, ...)
    FixedSizeArrays.wrap(FixedSizeArray, parent::DenseVector)

Return a [`FixedSizeArray`](@ref) of the given `size` which uses `parent` as its memory backend,
without copying.
`length(parent)` must be equal to `prod(size)`, otherwise a `DimensionMismatch` is thrown.
Without `size`, the result is a vector of length `length(parent)`.

The first argument can also be a more specific `FixedSizeArray` type, e.g.
`FixedSizeArray{T,N,Mem}`: the element type, the number of dimensions and the memory backend
are then checked against `parent` and `size`.

The result shares its elements with `parent`, and `parent(FixedSizeArrays.wrap(FixedSizeArray, p, size)) === p`.
Changing the length of `parent` (e.g. with `resize!` when `parent` is a `Vector`) while the
result is in use is not allowed: wrap `parent` again after resizing it instead.

This is similar to the (non-public) `Base.wrap(Array, mem, dims)` in Julia v1.11+,
except that `parent` must have exactly the length of the result.

---

    FixedSizeArrays.wrap(FixedSizeArray, ref::GenericMemoryRef, size::NTuple{N,Integer})
    FixedSizeArrays.wrap(FixedSizeArray, ref::GenericMemoryRef, size1::Integer, size2::Integer, ...)
    FixedSizeArrays.wrap(FixedSizeArray{T,N,MemoryRef{T}} where {T,N}, parent::Vector, size...)

Julia v1.11 and later: return a [`FixedSizeArray`](@ref) of the given `size` with a
`GenericMemoryRef` memory backend, whose elements are the `prod(size)` elements of the memory
of `ref` starting at `ref`. The memory must have at least that many elements from `ref` on,
otherwise a `DimensionMismatch` is thrown. `parent` of the result is `ref`.

With a `Vector` and a `FixedSizeArray` type with a `GenericMemoryRef` backend as the first
argument, the current reference to the first element of the `Vector` is wrapped. As for other
parents, the length of the `Vector` must be equal to `prod(size)`.

Unlike a `Vector` backend, the memory reference is an immutable field of the result, so the
pointer to the data does not need to be reloaded after writes to other arrays (e.g. in loops).
The result keeps the memory alive and stays valid when the `Vector` is resized, but it then
keeps referring to the memory the `Vector` had when it was wrapped: after `resize!` (or
`push!`, `popfirst!`, `sizehint!`, ...) of the `Vector` the result may no longer share its
elements with the `Vector`, and elements dropped by the `Vector` may be overwritten or
unassigned. Wrap the `Vector` again after resizing it.

```julia-repl
julia> v = collect(1:8); popfirst!(v);

julia> a = FixedSizeArrays.wrap(FixedSizeArray, v.ref, (2, 3));

julia> a, parent(a) === v.ref
([2 4 6; 3 5 7], true)

julia> b = FixedSizeArrays.wrap(FixedSizeArray{T,N,MemoryRef{T}} where {T,N}, v);

julia> typeof(b)
FixedSizeVector{Int64, MemoryRef{Int64}}
```

# Examples

```jldoctest
julia> using FixedSizeArrays

julia> v = collect(1:6);

julia> a = FixedSizeArrays.wrap(FixedSizeArray, v, (2, 3));

julia> size(a), parent(a) === v
((2, 3), true)

julia> a[2, 3] = 60;  # writes to `v`

julia> v[6]
60

julia> resize!(v, 8);  # invalidates `a`

julia> b = FixedSizeArrays.wrap(FixedSizeArray, v, 2, 4);  # wrap again after resizing

julia> size(b), b[2, 3]
((2, 4), 60)
```
"""
function wrap(::Type{FSA}, parent::DenseVector, size::Tuple{Vararg{Integer}}) where {FSA <: FixedSizeArray}
    size = check_ndims(FSA, size)
    s = map(Int, size)
    len = checked_dims(s)
    if length(parent) != len
        throw(DimensionMismatch(lazy"length of the parent ($(length(parent))) does not match the product of the dimensions $s"))
    end
    a = new_fixed_size_array(wrap_backend(FSA, parent), s)
    if !(a isa FSA)
        throw(ArgumentError(lazy"a `FixedSizeArray` wrapping a `$(typeof(parent))` with $(length(s)) dimensions is not a `$FSA`"))
    end
    a
end
wrap_backend(::Type{<:FixedSizeArray}, parent) = parent
function wrap(::Type{FSA}, parent::DenseVector, size1::Integer, size::Vararg{Integer}) where {FSA <: FixedSizeArray}
    wrap(FSA, parent, (size1, size...))
end
function wrap(::Type{FSA}, parent::DenseVector) where {FSA <: FixedSizeArray}
    wrap(FSA, parent, (length(parent),))
end

# `iterate`: the `AbstractArray` fallback doesn't perform well, so add our own methods (copied from `Base.Array`)

Base.iterate(A::FixedSizeArray, i=1) = (@inline; (i - 1)%UInt < length(A)%UInt ? (@inbounds A[i], i + 1) : nothing)

# destructuring interface: `Base.rest`: add a method of our own, improving performance

function Base.rest(a::FixedSizeArray, i::Int = 1)
    a[i:end]
end

"""
    FixedSizeArrayDefault{T,N}(undef, size::NTuple{N,Int})
    FixedSizeArrayDefault{T,N}(undef, size1::Int, size2::Int, ...)
    FixedSizeArrayDefault{T,N}(array::AbstractArray)

Construct a [`FixedSizeArray`](@ref) with element type `T`, number of dimensions `N`, and the default memory backend (`Vector{T}` on Julia v1.10, `Memory{T}` on Julia v1.11+).
"""
const FixedSizeArrayDefault = FixedSizeArray{T, N, default_underlying_storage_type{T}} where {T, N}
"""
    FixedSizeVectorDefault{T}(undef, size::NTuple{1,Int})
    FixedSizeVectorDefault{T}(undef, size1::Int)
    FixedSizeVectorDefault{T}(array::AbstractVector)

Construct a [`FixedSizeVector`](@ref) with element type `T`, and the default memory backend (`Vector{T}` on Julia v1.10, `Memory{T}` on Julia v1.11+).
"""
const FixedSizeVectorDefault = FixedSizeArrayDefault{T, 1} where {T}
"""
    FixedSizeMatrixDefault{T}(undef, size::NTuple{2,Int})
    FixedSizeMatrixDefault{T}(undef, size1::Int, size2::Int)
    FixedSizeMatrixDefault{T}(array::AbstractMatrix)

Construct a [`FixedSizeMatrix`](@ref) with element type `T`, and the default memory backend (`Vector{T}` on Julia v1.10, `Memory{T}` on Julia v1.11+).
"""
const FixedSizeMatrixDefault = FixedSizeArrayDefault{T, 2} where {T}
