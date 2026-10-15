"""
    static_dims(a::StaticSizeArray)::Tuple

Return the dimensions of `a`, a tuple of `StaticInt`s and `Int`s.
"""
static_dims(a::StaticSizeArray) = getfield(a, :dims)

# `Integer`s as `Int`s, `Colon`s (and `StaticInt`s, in the Static.jl extension) as they are
to_static_dim(d::Colon) = d
to_static_dim(d::Integer) = Int(d)::Int

Base.size(a::StaticSizeArray) = map(dim_value, static_dims(a))

# a method of our own, so that e.g. `size(a, 1)` is a compile-time constant for a static
# dimension (from `known_dims(D)`, a tuple of constants, as indexing the dimensions
# themselves with a `d` known only at runtime would box them)
Base.@constprop :aggressive function Base.size(a::StaticSizeArray{<:Any, N, D}, d::Integer) where {N, D}
    if d < 1
        throw(ArgumentError("dimension out of range"))
    elseif d ≤ N
        k = known_dims(D)[d]
        (k === nothing) ? size(a)[d] : k
    else
        1
    end
end

parent_type(::Type{<:StaticSizeArray{<:Any,<:Any,<:Any,Mem}}) where {Mem} = Mem

# `wrap`

function check_static_size_type(::Type{SA}, a::StaticSizeArray) where {SA <: StaticSizeArray}
    if !(a isa SA)
        throw(ArgumentError(lazy"a `$(typeof(a))` is not a `$SA`"))
    end
    a
end

# Replace a `Colon` in `dims` by the size computed from `len`.
function resolve_colon(dims::Tuple, len::Int)
    any(d -> d isa Colon, dims) || return dims
    others = filter(d -> !(d isa Colon), dims)
    if length(others) != length(dims) - 1
        throw(ArgumentError(lazy"at most one dimension can be a `Colon`, got $dims"))
    end
    m = checked_dims(map(dim_value, others))
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
    FixedSizeArrays.wrap(StaticSizeArray, parent::DenseVector, dims::Tuple)

Return a [`StaticSizeArray`](@ref) with the dimensions `dims`, a tuple of `StaticInt`s and
`Integer`s, which uses `parent` as its memory backend, without copying.

`length(parent)` must be equal to the product of `dims`, otherwise a `DimensionMismatch` is
thrown. At most one dimension can be given as `:`, it is then computed from
`length(parent)`. The length of `parent` must not be changed while the result is in use.

The first argument can also be a more specific type, e.g. `StaticSizeArray{T,N,D,Mem}`,
which is then checked against the result.

# Examples

```jldoctest
julia> using FixedSizeArrays: FixedSizeArrays, StaticSizeArray

julia> using Static: static

julia> v = collect(1.0:12.0);

julia> a = FixedSizeArrays.wrap(StaticSizeArray, v, (static(2), 3, static(2)));

julia> size(a), parent(a) === v
((2, 3, 2), true)

julia> size(FixedSizeArrays.wrap(StaticSizeArray, resize!(v, 18), (static(2), static(3), :)))  # wrap again after resizing
(2, 3, 3)
```
"""
function wrap(::Type{SA}, parent::DenseVector, dims::Tuple) where {SA <: StaticSizeArray}
    d = resolve_colon(map(to_static_dim, dims), length(parent))
    len = checked_dims(map(dim_value, d))
    if length(parent) != len
        throw(DimensionMismatch(lazy"length of the parent ($(length(parent))) does not match the product of the dimensions $d"))
    end
    check_static_size_type(SA, new_static_size_array(parent, d))
end

# constructors

# The dimensions type `D` and the memory backend of `SA`, `nothing` and the default memory
# backend respectively if they are not given, e.g. for `StaticSizeArray{T}`.
dims_type(::Type{<:StaticSizeArray{<:Any, <:Any, D}}) where {D <: Tuple} = D
dims_type(::Type{<:StaticSizeArray}) = nothing
parent_type_with_default(::Type{<:StaticSizeArray{T, <:Any, <:Any, Mem}}) where {T, Mem <: DenseVector{T}} = Mem
parent_type_with_default(::Type{<:StaticSizeArray{T}}) where {T} = default_underlying_storage_type{T}

# Convert the dimensions `dims` to the type `D`: the entries which are static in `D` must
# be equal to the static dimensions, the others are converted to `Int`s.
convert_dims(::Nothing, dims::Tuple) = map(to_static_dim, dims)
function convert_dims(::Type{D}, dims::Tuple) where {D <: Tuple}
    if fieldcount(D) != length(dims)
        throw(DimensionMismatch(lazy"the dimensions $dims do not match the dimensions type $D"))
    end
    map(convert_dim, fieldtypes(D), dims)
end

# (static dimensions: in the Static.jl extension)
convert_dim(::Type{Int}, d) = Int(dim_value(to_static_dim(d)))::Int

function (::Type{SA})(::UndefInitializer, dims::Tuple) where {T, SA <: StaticSizeArray{T}}
    d = convert_dims(dims_type(SA), dims)
    mem = parent_type_with_default(SA)(undef, checked_dims(map(dim_value, d)))
    check_static_size_type(SA, new_static_size_array(mem, d))
end

function StaticSizeArray(a::FixedSizeArray, dims::Tuple)
    if size(a) != map(dim_value ∘ to_static_dim, dims)
        throw(DimensionMismatch(lazy"the size of the array $(size(a)) does not match the dimensions $dims"))
    end
    new_static_size_array(parent(a), map(to_static_dim, dims))
end

# indexing (with one linear index: shared with `FixedSizeArray`s)

# The column-major linear index of the (in-bounds) indices `inds` for the dimensions `dims`,
# as a sum of strided terms, e.g. `1 + (i - 1) + n1 * (j - 1) + (n1 * n2) * (k - 1)` for the
# indices `(i, j, k)`, where the strides are constants as long as the dimensions before
# them are static.  (Unlike the Horner form `i + n1 * ((j - 1) + n2 * (k - 1))`, the terms
# of indices that don't change in a loop can be hoisted out of it if the dimensions are
# dynamic.)
function static_size_linear_index(dims::Tuple{Vararg{Any, N}}, inds::NTuple{N, Int}) where {N}
    terms = map((s, i) -> s * (i - 1), Base.size_to_strides(1, map(dim_value, dims)...), inds)
    +(1, terms...)
end

# The linear index of the (in-bounds) indices `inds`: with the static dimensions if there is
# one index per dimension, with the generic method otherwise (e.g. trailing ones).
function static_size_linear_index(a::StaticSizeArray{<:Any, N}, inds::Tuple{Vararg{Int, M}}) where {N, M}
    if M == N
        static_size_linear_index(static_dims(a), inds)
    else
        @inbounds LinearIndices(a)[inds...]
    end
end

Base.@propagate_inbounds function Base.getindex(a::StaticSizeArray, i1::Int, i2::Int, is::Vararg{Int, M}) where {M}
    inds = (i1, i2, is...)
    @boundscheck checkbounds(a, inds...)
    @inbounds parent(a)[static_size_linear_index(a, inds)]
end
Base.@propagate_inbounds function Base.getindex(a::StaticSizeArray{<:Any, N}, i::CartesianIndex{N}) where {N}
    a[Tuple(i)...]
end

Base.@propagate_inbounds @assume_noub_if_noinbounds function Base.setindex!(a::StaticSizeArray, x, i1::Int, i2::Int, is::Vararg{Int, M}) where {M}
    inds = (i1, i2, is...)
    @boundscheck checkbounds(a, inds...)
    @inbounds parent(a)[static_size_linear_index(a, inds)] = x
    a
end
Base.@propagate_inbounds function Base.setindex!(a::StaticSizeArray{<:Any, N}, x, i::CartesianIndex{N}) where {N}
    a[Tuple(i)...] = x
end

# strided arrays interface (`elsize`, `cconvert` and `unsafe_convert`: shared with
# `FixedSizeArray`s)

Base.strides(a::StaticSizeArray) = Base.size_to_strides(1, size(a)...)

# `similar`: a `StaticSizeArray` with the same dimensions without a size, a `FixedSizeArray`
# with a size (deciding by the value of the size would not be type-stable, e.g. for range
# indexing and reductions); `similar(typeof(a), size)` keeps the static dimensions

function Base.similar(a::StaticSizeArray, ::Type{E}) where {E}
    new_static_size_array(similar(parent(a), E, (length(a),)), static_dims(a))
end

function Base.similar(a::StaticSizeArray, ::Type{E}, size::NTuple{N, Int}) where {E, N}
    new_fixed_size_array(similar(parent(a), E, (checked_dims(size),)), size)
end

# `copy` (`fill!`, `reshape` and `copyto!`: shared with `FixedSizeArray`s)

function Base.copy(a::StaticSizeArray)
    new_static_size_array(copy(parent(a)), static_dims(a))
end

# broadcasting: as for `FixedSizeArray`s, the result is a `FixedSizeArray`

function Base.BroadcastStyle(::Union{Type{<:(StaticSizeArray{T, N, D, Mem} where {T, D})},
                                     Type{<:(SubArray{T, N, <:StaticSizeArray{T, N, <:Any, Mem}} where {T})}}) where {N, Mem}
    fixed_size_array_broadcast_style(Val(N), Mem)
end
