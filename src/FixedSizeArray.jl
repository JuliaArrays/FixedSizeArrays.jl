export
    BoundsErrorLight,
    FixedSizeArray, FixedSizeVector, FixedSizeMatrix,
    FixedSizeArrayDefault, FixedSizeVectorDefault, FixedSizeMatrixDefault

# Dimensions of a `FixedSizeArray`: each one is either dynamic, an `Int` stored in the
# array, or static, a compile-time constant which takes no space, e.g.
# `Tuple{StaticInt{4}, Int}` for `(4, n)`. Static dimensions are `StaticInt`s from
# Static.jl, which is a weak dependency: the methods for them are in the package extension
# `FixedSizeArraysStaticExt`, while `src` only uses these functions.

# The value of a dimension.
#
# An `Int` dimension is passed through an intrinsic, so that `size` loads all dimensions at
# once. With a plain `getfield(a, :dims)`, Julia's codegen loads each dimension of a
# by-reference `FixedSizeArray` argument lazily, at its use. In a short-circuiting check
# such as `axes(a) == expected || throw(...)` the loads of dimensions 2:N then only happen
# in the conditional blocks of the `&&` chain, so they do not dominate later uses in loops.
# LLVM's LICM later hoists the loads of the loop out of conditionally executed inner loops
# and drops their `!tbaa`/`!invariant.load` metadata, after which GVN can neither CSE them
# with the check nor use the checked values, so the dimensions stay dynamic in the loop.
# Passing each dimension through an intrinsic forces the tuple to be materialized in `size`
# (as for `Array` and `PtrArray`) and is a no-op otherwise.
dim_value(d::Int) = Core.Intrinsics.bitcast(Int, d)

# The value of a dimension of type `D` if it is static, `nothing` otherwise.
known_dim(::Type{Int}) = nothing
known_dims(::Type{D}) where {D <: Tuple} = map(known_dim, fieldtypes(D))

# (`D` has no bound in the type, so that e.g. `FixedSizeArrayDims{T, N, D, Mem} where {D}` is
# a `FixedSizeArrayDims`; the inner constructor only accepts `N` dimensions. The docstring is
# the one of `FixedSizeArray`, the alias with all dimensions dynamic.)
"""
    FixedSizeArrayDims{T,N,D,Mem<:DenseVector{T}} <: DenseArray{T,N}

The underlying type of [`FixedSizeArray`](@ref), with the dimensions type `D`, a tuple type
of `Int`s and `StaticInt`s from [Static.jl](https://github.com/SciML/Static.jl):
`FixedSizeArray{T,N,Mem}` is `FixedSizeArrayDims{T,N,NTuple{N,Int},Mem}`, with all
dimensions dynamic. Use `FixedSizeArrayDims` to construct a `FixedSizeArray` with static
dimensions, see [`FixedSizeArray`](@ref) and [`FixedSizeArrays.wrap`](@ref):

    FixedSizeArrayDims{T}(undef, dims::Tuple)
    FixedSizeArrayDims{T,N,D,Mem}(undef, dims::Tuple)
    FixedSizeArrayDims{T,N,D,Mem}(array::AbstractArray)
    FixedSizeArrays.wrap(FixedSizeArrayDims, parent::DenseVector, dims::Tuple)

!!! note

    `FixedSizeArrayDims` is public, but not exported. Use it as
    `FixedSizeArrays.FixedSizeArrayDims`, or import it with
    `using FixedSizeArrays: FixedSizeArrayDims`.
"""
struct FixedSizeArrayDims{T,N,D,Mem<:DenseVector{T}} <: DenseArray{T,N}
    mem::Mem
    dims::D
    global function new_fixed_size_array(mem::DenseVector{T}, dims::Tuple{Vararg{Any, N}}) where {T,N}
        new{T,N,typeof(dims),typeof(mem)}(mem, dims)
    end
end

"""
    FixedSizeArray{T,N,Mem<:DenseVector{T}}(undef, size::NTuple{N,Int})
    FixedSizeArray{T,N,Mem<:DenseVector{T}}(undef, size1::Int, size2::Int, ...)
    FixedSizeArray{T,N,Mem<:DenseVector{T}}(array::AbstractArray)

Construct a fixed-size `FixedSizeArray` with element type `T`, number of dimensions `N`, and memory backend `Mem`.

Convenient aliases provided:

* [`FixedSizeVector{T,Mem}`](@ref): for 1-dimensional fixed-size arrays
* [`FixedSizeMatrix{T,Mem}`](@ref): for 2-dimensional fixed-size arrays
* [`FixedSizeArrayDefault{T,N}`](@ref): for `N`-dimensional fixed-size arrays, with default memory backend
* [`FixedSizeVectorDefault{T}`](@ref): for 1-dimensional fixed-size arrays, with default memory backend
* [`FixedSizeMatrixDefault{T}`](@ref): for 2-dimensional fixed-size arrays, with default memory backend

# Static dimensions

`FixedSizeArray{T,N,Mem}` is an alias for [`FixedSizeArrayDims{T,N,NTuple{N,Int},Mem}`](@ref FixedSizeArrayDims):
the dimensions of a `FixedSizeArrayDims{T,N,D,Mem}` are a tuple `dims::D` of `Int`s and
`StaticInt`s from [Static.jl](https://github.com/SciML/Static.jl). Each dimension can be
dynamic, i.e. an `Int` stored in the array, or static, i.e. a compile-time constant which
is part of the type `D` and takes no space in the array. E.g. with
`dims = (static(5), static(4), n)`, `D = Tuple{StaticInt{5}, StaticInt{4}, Int}`, the array
has size `(5, 4, n)`, and `n` can be different for different arrays of the same type.
Static.jl is a weak dependency: load it (e.g. with `using Static: static`) to use static
dimensions.

The static dimensions allow the compiler to compute the linear indices, the strides, and
the sizes in these dimensions at compile time, which is useful e.g. for arrays where some
dimensions are fixed (variables, nodes, ...) while others (elements of a mesh, particles,
...) change at runtime. An array with static dimensions is constructed with the
[`FixedSizeArrayDims`](@ref) constructors, with `collect_as`, or with
[`FixedSizeArrays.wrap`](@ref) around an existing memory backend.

The dimensions given to a constructor are converted to the dimensions type of the array, if
its type gives it: e.g. `FixedSizeArray{T}(undef, (static(2), 3))` has two dynamic
dimensions, while `FixedSizeArrayDims{T}(undef, (static(2), 3))` has a static one.

# Examples

```jldoctest
julia> using FixedSizeArrays

julia> using FixedSizeArrays: FixedSizeArrayDims

julia> using Static: static, StaticInt

julia> a = FixedSizeArray{Float64}(undef, (2, 3));

julia> typeof(a) === FixedSizeArray{Float64, 2, Memory{Float64}} === FixedSizeArrayDims{Float64, 2, Tuple{Int, Int}, Memory{Float64}}
true

julia> b = FixedSizeArrayDims{Int}(undef, (static(2), static(3), 2));  # two static dimensions

julia> typeof(b)
FixedSizeArrayDims{Int64, 3, Tuple{StaticInt{2}, StaticInt{3}, Int64}, Memory{Int64}}

julia> size(b), size(b, 1), strides(b)
((2, 3, 2), 2, (1, 2, 6))
```
"""
const FixedSizeArray{T,N,Mem<:DenseVector{T}} = FixedSizeArrayDims{T,N,NTuple{N,Int},Mem}

function Base.propertynames(
    # the `unused` is here because of https://github.com/JuliaLang/julia/issues/44428
    (@nospecialize unused::FixedSizeArrayDims),
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
    special_storage_types = (Vector, optional_memory...)
    Union{
        # 0 fixed parameters
        Type{FixedSizeArray},
        # 1 fixed parameter
        Type{FixedSizeArray{T}} where {T},
        Type{FixedSizeArray{T, N} where {T}} where {N},
        # 2 fixed parameters
        Type{FixedSizeArray{T, N}} where {T, N},
        Type{FixedSizeArray{T, N, Mem} where {N}} where {T, Mem <: DenseVector{T}},
        # 3 fixed parameters
        Type{FixedSizeArray{T, N, Mem}} where {T, N, Mem <: DenseVector{T}},
        # special cases depending on the underlying storage type
        map(f, special_storage_types)...,
    }
end

# (only the `FixedSizeArray` constructors are restricted)
check_constructor_is_allowed(::Type{T}) where {T <: FixedSizeArrayDims} = T
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

"""
    static_dims(a::FixedSizeArray)::Tuple

Return the dimensions of `a`, a tuple of `StaticInt`s and `Int`s.
"""
static_dims(a::FixedSizeArrayDims) = getfield(a, :dims)

# `Integer`s as `Int`s, `StaticInt`s (in the Static.jl extension) as they are
to_static_dim(d::Integer) = Int(d)::Int
to_static_dim(::Colon) = throw(ArgumentError("a `Colon` dimension is only supported by `FixedSizeArrays.wrap`"))

Base.size(a::FixedSizeArrayDims) = map(dim_value, static_dims(a))

# a method of our own, so that e.g. `size(a, 1)` is a compile-time constant for a static
# dimension (from `known_dims(D)`, a tuple of constants, as indexing the dimensions
# themselves with a `d` known only at runtime would box them)
Base.@constprop :aggressive function Base.size(a::FixedSizeArrayDims{<:Any, N, D}, d::Integer) where {N, D}
    if d < 1
        throw(ArgumentError("dimension out of range"))
    elseif d ≤ N
        k = known_dims(D)[d]
        (k === nothing) ? size(a)[d] : k
    else
        1
    end
end

parent_type(::Type{<:FixedSizeArrayDims{<:Any,<:Any,<:Any,Mem}}) where {Mem} = Mem

# The memory backend of `FSA`, the default memory backend if it is not given.
function parent_type_with_default(::Type{<:(FixedSizeArrayDims{T, N, D, Mem} where {N, D})}) where {T, Mem <: DenseVector{T}}
    Mem
end
function parent_type_with_default(::Type{<:FixedSizeArrayDims{T}}) where {T}
    default_underlying_storage_type{T}
end
function parent_type_with_default(::Type{<:FixedSizeArrayDims})
    default_underlying_storage_type
end
for Mem ∈ (Vector, optional_memory...)
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

# The dimensions type `D` of `FSA`, if it is given, `nothing` otherwise (e.g. for
# `FixedSizeArrayDims{T}`). For a `FixedSizeArray` without `N` it is `Tuple{Vararg{Int}}`:
# all dimensions are `Int`s.
dims_type(::Type{<:FixedSizeArrayDims{<:Any, <:Any, D}}) where {D <: Tuple} = D
dims_type(::Type{<:FixedSizeArray}) = Tuple{Vararg{Int}}
dims_type(::Type{<:FixedSizeArrayDims}) = nothing

# Check that the number of dimensions matches `N`, if `FSA` gives it.
function check_ndims(::Type{<:(FixedSizeArrayDims{T, N} where {T})}, dims::Tuple) where {N}
    if dims isa Tuple{Vararg{Any, N}}
        dims
    else
        throw(DimensionMismatch("mismatch between dimension count in type and the length of the size tuple"))
    end
end
function check_ndims(::Type{<:FixedSizeArrayDims}, dims::Tuple)
    dims
end

# Convert the dimensions `dims` to the type `D`: the entries which are static in `D` must
# be equal to the static dimensions, the others are converted to `Int`s. Without `D`, the
# dimensions are taken as given: `Integer`s as `Int`s, `StaticInt`s as they are.
convert_dims(::Nothing, dims::Tuple) = map(to_static_dim, dims)
convert_dims(::Type{Tuple{Vararg{Int}}}, dims::Tuple) = map(Base.Fix1(convert_dim, Int), dims)
function convert_dims(::Type{D}, dims::Tuple) where {D <: Tuple}
    if fieldcount(D) != length(dims)
        throw(DimensionMismatch(lazy"the dimensions $dims do not match the dimensions type $D"))
    end
    map(convert_dim, fieldtypes(D), dims)
end

# (static dimensions: in the Static.jl extension)
convert_dim(::Type{Int}, d) = Int(dim_value(to_static_dim(d)))::Int

# The dimensions `dims` of an `FSA`, converted to its dimensions type.
function to_dims(::Type{FSA}, dims::Tuple) where {FSA <: FixedSizeArrayDims}
    convert_dims(dims_type(FSA), check_ndims(FSA, dims))
end

function check_array_type(::Type{FSA}, a::FixedSizeArrayDims) where {FSA <: FixedSizeArrayDims}
    if !(a isa FSA)
        throw(ArgumentError(lazy"a `$(typeof(a))` is not a `$FSA`"))
    end
    a
end

# `undef` constructors

function undef_constructor(::Type{FSA}, dims::Tuple) where {T, FSA <: FixedSizeArrayDims{T}}
    d = to_dims(FSA, dims)
    mem = parent_type_with_default(FSA)(undef, checked_dims(map(dim_value, d)))
    check_array_type(FSA, new_fixed_size_array(mem, d))
end

function (::Type{FSA})(::UndefInitializer, dims::Tuple) where {T, FSA <: FixedSizeArrayDims{T}}
    undef_constructor(FSA, dims)
end
function (::Type{FSA})(::UndefInitializer, dims::Vararg{Integer}) where {T, FSA <: FixedSizeArrayDims{T}}
    undef_constructor(FSA, dims)
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

# indexing

Base.IndexStyle(::Type{<:FixedSizeArrayDims}) = IndexLinear()
Base.@propagate_inbounds function Base.getindex(A::FixedSizeArrayDims, i::Int)
    @boundscheck checkbounds(A, i)
    @inbounds parent(A)[i]
end
Base.@propagate_inbounds @assume_noub_if_noinbounds function Base.setindex!(A::FixedSizeArrayDims, x, i::Int)
    @boundscheck checkbounds(A, i)
    @inbounds parent(A)[i] = x
    return A
end

function Base.checkbounds(A::FixedSizeArrayDims, is...)
    checkbounds_lightboundserror(A, is...)
end

function Base.isassigned(a::FixedSizeArrayDims, i::Int)
    @boundscheck checkbounds(Bool, a, i) || return false
    @inbounds isassigned(parent(a), i)
end

# Cartesian indexing: the column-major linear index of the (in-bounds) indices `inds` for
# the dimensions `dims`, as a sum of strided terms, e.g.
# `1 + (i - 1) + n1 * (j - 1) + (n1 * n2) * (k - 1)` for the indices `(i, j, k)`, where the
# strides are constants as long as the dimensions before them are static. This is the
# computation of `Base` (`_sub2ind`), but with the strides from the static dimensions.
function static_size_linear_index(dims::Tuple{Vararg{Any, N}}, inds::NTuple{N, Int}) where {N}
    terms = map((s, i) -> s * (i - 1), Base.size_to_strides(1, map(dim_value, dims)...), inds)
    +(1, terms...)
end

# The linear index of the (in-bounds) indices `inds`: with the dimensions if there is one
# index per dimension, with the generic method otherwise (e.g. trailing ones).
function static_size_linear_index(a::FixedSizeArrayDims{<:Any, N}, inds::Tuple{Vararg{Int, M}}) where {N, M}
    if M == N
        static_size_linear_index(static_dims(a), inds)
    else
        @inbounds LinearIndices(a)[inds...]
    end
end

Base.@propagate_inbounds function Base.getindex(a::FixedSizeArrayDims, i1::Int, i2::Int, is::Vararg{Int, M}) where {M}
    inds = (i1, i2, is...)
    @boundscheck checkbounds(a, inds...)
    @inbounds parent(a)[static_size_linear_index(a, inds)]
end
Base.@propagate_inbounds function Base.getindex(a::FixedSizeArrayDims{<:Any, N}, i::CartesianIndex{N}) where {N}
    a[Tuple(i)...]
end

Base.@propagate_inbounds @assume_noub_if_noinbounds function Base.setindex!(a::FixedSizeArrayDims, x, i1::Int, i2::Int, is::Vararg{Int, M}) where {M}
    inds = (i1, i2, is...)
    @boundscheck checkbounds(a, inds...)
    @inbounds parent(a)[static_size_linear_index(a, inds)] = x
    a
end
Base.@propagate_inbounds function Base.setindex!(a::FixedSizeArrayDims{<:Any, N}, x, i::CartesianIndex{N}) where {N}
    a[Tuple(i)...] = x
end

# `similar`: an array with the same dimensions (also the static ones) without a size, a
# `FixedSizeArray` (all dimensions dynamic) with a size (deciding by the value of the size
# would not be type-stable, e.g. for range indexing and reductions);
# `similar(typeof(a), size)` keeps the static dimensions

function Base.similar(a::FixedSizeArrayDims, ::Type{E}) where {E}
    new_fixed_size_array(similar(parent(a), E, (length(a),)), static_dims(a))
end

function Base.similar(a::FixedSizeArrayDims, ::Type{E}, size::NTuple{N, Int}) where {E, N}
    new_fixed_size_array(similar(parent(a), E, (checked_dims(size),)), size)
end

# strided arrays interface (`elsize`, `cconvert` and `unsafe_convert`: below)

Base.strides(a::FixedSizeArrayDims) = Base.size_to_strides(1, size(a)...)

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

# broadcasting: the result of out-of-place broadcasting is a `FixedSizeArray`, also for
# arrays with static dimensions

struct FixedSizeArrayBroadcastStyle{N, Mem <: DenseVector} <: Broadcast.AbstractArrayStyle{N} end

function (::Type{<:(FixedSizeArrayBroadcastStyle{N, Mem} where {N})})(::Val{M}) where {Mem, M}
    m = check_count_value(M)
    FixedSizeArrayBroadcastStyle{m, Mem}()  # `FixedSizeArray` supports any dimensionality
end

function Base.BroadcastStyle(::Union{Type{<:(FixedSizeArrayDims{T, N, D, Mem} where {T, D})},
                                     Type{<:(SubArray{T, N, <:FixedSizeArrayDims{T, N, <:Any, Mem}} where {T})}}) where {N, Mem}
    fixed_size_array_broadcast_style(Val(N), Mem)
end

function fixed_size_array_broadcast_style(::Val{N}, ::Type{Mem}) where {N, Mem}
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

"""
    parent(f::FixedSizeArray)

Return the underlying parent object of the fixed-size array `f`, a `DenseVector` of length
`length(f)`. Use `parentindices(f)` to get the valid indices of the parent of `f`.
"""
Base.parent(f::FixedSizeArrayDims) = getfield(f, :mem)

Base.parentindices(f::FixedSizeArrayDims) = axes(parent(f))

underlying_storage(m) = m
underlying_storage(f::FixedSizeArrayDims) = parent(f)

axes_are_one_based(axes) = all(isone ∘ first, axes)

# converting constructors for copying other array types

function converting_constructor(::Type{FSA}, src::AbstractArray) where {FSA <: FixedSizeArrayDims}
    axs = axes(src)
    if !axes_are_one_based(axs)
        throw(DimensionMismatch("source array has a non-one-based indexing axis"))
    end
    collect_as_haseltype(FSA, src)
end

function (::Type{FSA})(src::AbstractArray) where {N, FSA <: FixedSizeArrayDims{<:Any, N}}
    if N::Int != ndims(src)
        throw(DimensionMismatch("dimensionality mismatch"))
    end
    converting_constructor(FSA, src)
end

function (::Type{FSA})(src::AbstractArray) where {FSA <: FixedSizeArrayDims}
    converting_constructor(FSA, src)
end

# `copy`: avoid the `similar` and `copyto!` in the generic fallback

function Base.copy(a::FixedSizeArrayDims)
    new_fixed_size_array(copy(parent(a)), static_dims(a))
end

# `getindex` with a `Colon` for the index: better effects than with the generic fallback

function Base.getindex(a::FixedSizeArrayDims, ::Colon)
    c = copy(a)
    if a isa AbstractVector
        c
    else
        reshape(c, :)
    end
end

# `copyto!`

Base.@propagate_inbounds function copyto5!(dst, doff, src, soff, n)
    copyto!(underlying_storage(dst), doff, underlying_storage(src), soff, n)
    dst
end

Base.@propagate_inbounds function copyto2!(dst, src)
    copyto!(underlying_storage(dst), underlying_storage(src))
    dst
end

Base.@propagate_inbounds Base.copyto!(dst::FixedSizeArrayDims, doff::Integer, src::FixedSizeArrayDims, soff::Integer, n::Integer) = copyto5!(dst, doff, src, soff, n)
Base.@propagate_inbounds Base.copyto!(dst::FixedSizeArrayDims, src::FixedSizeArrayDims) = copyto2!(dst, src)

for A ∈ (Array, optional_generic_memory...)  # Add more? Too bad we have to hardcode to avoid ambiguity.
    @eval begin
        Base.@propagate_inbounds Base.copyto!(dst::FixedSizeArrayDims, doff::Integer, src::$A,                 soff::Integer, n::Integer) = copyto5!(dst, doff, src, soff, n)
        Base.@propagate_inbounds Base.copyto!(dst::$A,                 doff::Integer, src::FixedSizeArrayDims, soff::Integer, n::Integer) = copyto5!(dst, doff, src, soff, n)

        Base.@propagate_inbounds Base.copyto!(dst::FixedSizeArrayDims, src::$A                ) = copyto2!(dst, src)
        Base.@propagate_inbounds Base.copyto!(dst::$A,                 src::FixedSizeArrayDims) = copyto2!(dst, src)
    end
end

# `fill!`: loop over the parent, as the generic fallback, iterating over `eachindex` of the
# array itself, does not elide the bounds checks on Julia v1.13 and later, so does not
# vectorize: <https://github.com/JuliaArrays/FixedSizeArrays.jl/issues/200>.
# Forwarding to `fill!` of the parent instead is slower for short arrays on these versions.

function Base.fill!(a::FixedSizeArrayDims{T}, x) where {T}
    xT = x isa T ? x : convert(T, x)::T
    mem = parent(a)
    for i ∈ eachindex(mem)
        @inbounds mem[i] = xT
    end
    a
end

# `copy!`: avoid the `resize!` branch of `copy!(::AbstractVector, ::AbstractVector)`

function Base.copy!(dst::FixedSizeArrayDims{<:Any, 1}, src::AbstractVector)
    axes(dst) == axes(src) || throw(ArgumentError(lazy"arrays must have the same axes for copy! (consider using `copyto!`), got $(axes(dst)) and $(axes(src))"))
    copyto!(dst, src)
end

# unsafe: the native address of the array's storage

Base.cconvert(::Type{<:Ptr}, a::FixedSizeArrayDims) = parent(a)

function Base.unsafe_convert(::Type{Ptr{T}}, a::FixedSizeArrayDims{T}) where {T}
    Base.unsafe_convert(Ptr{T}, parent(a))
end

# `elsize`: part of the strided arrays interface, used for `pointer`

Base.elsize(::Type{A}) where {A<:FixedSizeArrayDims} = Base.elsize(parent_type(A))

# `reshape`: specializing it to ensure it returns a `FixedSizeArray` (all dimensions dynamic)

function Base.reshape(a::FixedSizeArrayDims, size::(NTuple{N,Int} where {N}))
    len = checked_dims(size)
    if length(a) != len
        throw(DimensionMismatch("new shape not consistent with existing array length"))
    end
    new_fixed_size_array(parent(a), size)
end

# `wrap`: construct a `FixedSizeArray` around an existing memory backend, without copying

# Replace a `Colon` in `dims` by the size computed from `len`.
function resolve_colon(dims::Tuple, len::Int)
    any(d -> d isa Colon, dims) || return dims
    others = filter(d -> !(d isa Colon), dims)
    if length(others) != length(dims) - 1
        throw(ArgumentError(lazy"at most one dimension can be a `Colon`, got $dims"))
    end
    m = checked_dims(map(dim_value ∘ to_static_dim, others))
    if iszero(m)
        throw(ArgumentError(lazy"can not compute the size of the `Colon` dimension from the length of the parent for the dimensions $dims"))
    end
    (n, r) = divrem(len, m)
    if !iszero(r)
        throw(DimensionMismatch(lazy"length of the parent ($len) is not a multiple of the product of the dimensions $others"))
    end
    map(d -> d isa Colon ? n : d, dims)
end

"""
    FixedSizeArrays.wrap(FixedSizeArray, parent::DenseVector, size::NTuple{N,Integer})
    FixedSizeArrays.wrap(FixedSizeArray, parent::DenseVector, size1::Integer, size2::Integer, ...)
    FixedSizeArrays.wrap(FixedSizeArray, parent::DenseVector)
    FixedSizeArrays.wrap(FixedSizeArrayDims, parent::DenseVector, size::NTuple{N,Integer})
    FixedSizeArrays.wrap(FixedSizeArrayDims, parent::DenseVector, size1::Integer, size2::Integer, ...)
    FixedSizeArrays.wrap(FixedSizeArrayDims, parent::DenseVector)

Return a [`FixedSizeArray`](@ref) of the given `size` which uses `parent` as its memory backend,
without copying.
`length(parent)` must be equal to `prod(size)`, otherwise a `DimensionMismatch` is thrown.
At most one dimension can be given as `:`, it is then computed from `length(parent)`.
Without `size`, the result is a vector of length `length(parent)`.

With [`FixedSizeArrayDims`](@ref) as the first argument, the dimensions can also be
`StaticInt`s, for an array with static dimensions.

The first argument can also be a more specific type, e.g. `FixedSizeArray{T,N,Mem}` or
`FixedSizeArrayDims{T,N,D,Mem}`: the dimensions are converted to the dimensions type of the
array (to `Int`s for a `FixedSizeArray`), and the element type, the number of dimensions
and the memory backend are checked against `parent` and the dimensions.

The result shares its elements with `parent`, and `parent(FixedSizeArrays.wrap(FixedSizeArray, p, size)) === p`.
Changing the length of `parent` (e.g. with `resize!` when `parent` is a `Vector`) while the
result is in use is not allowed: wrap `parent` again after resizing it instead.

This is similar to the (non-public) `Base.wrap(Array, mem, dims)` in Julia v1.11+,
except that `parent` must have exactly the length of the result.

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

With static dimensions:

```jldoctest
julia> using FixedSizeArrays: FixedSizeArrays, FixedSizeArrayDims

julia> using Static: static

julia> v = collect(1.0:12.0);

julia> a = FixedSizeArrays.wrap(FixedSizeArrayDims, v, (static(2), 3, static(2)));

julia> size(a), parent(a) === v
((2, 3, 2), true)

julia> size(FixedSizeArrays.wrap(FixedSizeArrayDims, resize!(v, 18), (static(2), static(3), :)))  # wrap again after resizing
(2, 3, 3)
```
"""
function wrap(::Type{FSA}, parent::DenseVector, dims::Tuple) where {FSA <: FixedSizeArrayDims}
    d = convert_dims(dims_type(FSA), resolve_colon(check_ndims(FSA, dims), length(parent)))
    len = checked_dims(map(dim_value, d))
    if length(parent) != len
        throw(DimensionMismatch(lazy"length of the parent ($(length(parent))) does not match the product of the dimensions $d"))
    end
    check_array_type(FSA, new_fixed_size_array(parent, d))
end
function wrap(::Type{FSA}, parent::DenseVector, size1::Integer, size::Vararg{Integer}) where {FSA <: FixedSizeArrayDims}
    wrap(FSA, parent, (size1, size...))
end
function wrap(::Type{FSA}, parent::DenseVector) where {FSA <: FixedSizeArrayDims}
    wrap(FSA, parent, (length(parent),))
end

# `iterate`: the `AbstractArray` fallback doesn't perform well, so add our own methods (copied from `Base.Array`)

Base.iterate(A::FixedSizeArrayDims, i=1) = (@inline; (i - 1)%UInt < length(A)%UInt ? (@inbounds A[i], i + 1) : nothing)

# destructuring interface: `Base.rest`: add a method of our own, improving performance

function Base.rest(a::FixedSizeArrayDims, i::Int = 1)
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

# Printing of the types, as when `FixedSizeArray` was the `struct`: Julia prints a type
# through an alias if exactly one alias matches it, and `FixedSizeArray` is now an alias too,
# so it would make the other aliases ambiguous. Instead, a `FixedSizeArray` type is printed
# as one of the aliases below if exactly that one matches it, and as a `FixedSizeArray{T,N,Mem}`
# otherwise, independently of the Julia version. (This is a method of `show` rather than of
# `Base.make_typealias`, as a method of the latter invalidates much of the printing code.
# `make_wheres`, `show_typealias` & co. are internal functions of `Base`, used by `show`
# for types since Julia v1.6; if they disappear, only the printing changes.)
const printing_aliases = (
    :FixedSizeVector, :FixedSizeMatrix,
    :FixedSizeArrayDefault, :FixedSizeVectorDefault, :FixedSizeMatrixDefault,
)

# The parameters of the alias `alias` which give `x`, `nothing` if there are none (adapted
# from `Base.make_typealias`).
function alias_parameters(@nospecialize(x::Type), @nospecialize(alias::UnionAll))
    x <: alias || return nothing
    (_, env) = ccall(:jl_type_intersection_with_env, Any, (Any, Any), x, alias)::Core.SimpleVector
    applied = try
        alias{env...}
    catch ex
        ex isa TypeError || rethrow()
        return nothing
    end
    applied = Base.rewrap_unionall(applied, x)
    (!Base.has_free_typevars(applied) && applied === x) ? env : nothing
end

# The alias to print `x` with, and its parameters.
function printing_alias(@nospecialize(x::Type{<:FixedSizeArray}))
    match = nothing
    for name ∈ printing_aliases
        env = alias_parameters(x, getfield(@__MODULE__, name))
        env === nothing && continue
        match === nothing || return fixed_size_array_parameters(x)  # ambiguous
        match = (GlobalRef(@__MODULE__, name), env)
    end
    something(match, fixed_size_array_parameters(x))
end

# the parameters `T`, `N` and `Mem` of `FixedSizeArray` (the dimensions type is `NTuple{N,Int}`)
function fixed_size_array_parameters(@nospecialize(x::Type{<:FixedSizeArray}))
    body = Base.unwrap_unionall(x)::DataType
    (GlobalRef(@__MODULE__, :FixedSizeArray), Core.svec(body.parameters[1], body.parameters[2], body.parameters[4]))
end

if all(f -> isdefined(Base, f), (:makeproper, :make_wheres, :show_typealias, :show_wheres))
    function Base.show(io::IO, @nospecialize(x::Type{<:FixedSizeArray}))
        properx = Base.makeproper(io, x)
        # (not for `Union{}` and `Union`s of `FixedSizeArray` types, which are `FixedSizeArray`
        # types too, but have no parameters)
        if !(Base.unwrap_unionall(properx) isa DataType) || !(properx <: FixedSizeArray) || !(get(io, :compact, true)::Bool)
            return invoke(show, Tuple{IO, Type}, io, x)
        end
        # the aliases themselves are printed by their name, without a `where` clause
        for name ∈ (:FixedSizeArray, printing_aliases...)
            if x === getfield(@__MODULE__, name)
                return Base.show_typealias(io, GlobalRef(@__MODULE__, name), x, Core.svec(), TypeVar[])
            end
        end
        (name, env) = printing_alias(properx)
        wheres = Base.make_wheres(io, env, x)
        Base.show_typealias(io, name, x, env, wheres)
        Base.show_wheres(io, wheres)
    end
end
