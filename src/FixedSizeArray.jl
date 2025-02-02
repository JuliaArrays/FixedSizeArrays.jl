export
    FixedSizeArray, FixedSizeVector, FixedSizeMatrix,
    FixedSizeArrayDefault, FixedSizeVectorDefault, FixedSizeMatrixDefault

struct FixedSizeArray{T,N,Mem<:DenseVector{T}} <: DenseArray{T,N}
    mem::Mem
    size::NTuple{N,Int}
    global function new_fixed_size_array(mem::M, size::NTuple{N,Int}) where {T,N,M<:DenseVector{T}}
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

function FixedSizeArray{T,N,V}(::UndefInitializer, size::NTuple{N,Int}) where {T,N,V}
    new_fixed_size_array(V(undef, checked_dims(size))::V, size)
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
Base.@propagate_inbounds Base.getindex(A::FixedSizeArray, i::Int) = A.mem[i]
Base.@propagate_inbounds @assume_noub_if_noinbounds function Base.setindex!(A::FixedSizeArray, x, i::Int)
    @boundscheck checkbounds(A, i)
    @inbounds A.mem[i] = x
    return A
end

Base.size(a::FixedSizeArray) = a.size

function Base.similar(::T, ::Type{E}, size::NTuple{N,Int}) where {T<:FixedSizeArray,E,N}
    spec = TypeParametersElementTypeAndDimensionality()
    S = val_parameter(with_stripped_type_parameters(spec, T)){E, N}
    S(undef, size)
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

struct FixedSizeArrayBroadcastStyle{N, Mem <: DenseVector} <: Broadcast.AbstractArrayStyle{N} end

function (::Type{<:(FixedSizeArrayBroadcastStyle{N, Mem} where {N})})(::Val{M}) where {Mem, M}
    m = check_count_value(M)
    FixedSizeArrayBroadcastStyle{m, Mem}()  # `FixedSizeArray` supports any dimensionality
end

function Base.BroadcastStyle(::Type{<:(FixedSizeArray{T, N, Mem} where {T})}) where {N, Mem}
    n = check_count_value(N)
    spec = TypeParametersElementType()
    mem = val_parameter(with_stripped_type_parameters(spec, Mem))
    FixedSizeArrayBroadcastStyle{n, mem}()
end

function Base.similar(
    bc::Broadcast.Broadcasted{FixedSizeArrayBroadcastStyle{N, Mem}},
    ::Type{E},
) where {N,Mem,E}
    similar(FixedSizeArray{E, N, Mem{E}}, axes(bc))
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
    s = s::(Type{T} where {T>:t})
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

for A ∈ (Array, optional_generic_memory...)  # Add more? Too bad we have to hardcode to avoid ambiguity.
    @eval begin
        Base.@propagate_inbounds Base.copyto!(dst::FixedSizeArray, doff::Integer, src::$A,             soff::Integer, n::Integer) = copyto5!(dst, doff, src, soff, n)
        Base.@propagate_inbounds Base.copyto!(dst::$A,             doff::Integer, src::FixedSizeArray, soff::Integer, n::Integer) = copyto5!(dst, doff, src, soff, n)

        Base.@propagate_inbounds Base.copyto!(dst::FixedSizeArray, src::$A            ) = copyto2!(dst, src)
        Base.@propagate_inbounds Base.copyto!(dst::$A,             src::FixedSizeArray) = copyto2!(dst, src)
    end
end

# unsafe: the native address of the array's storage

Base.cconvert(::Type{<:Ptr}, a::FixedSizeArray) = a.mem

function Base.unsafe_convert(::Type{Ptr{T}}, a::FixedSizeArray{T}) where {T}
    Base.unsafe_convert(Ptr{T}, a.mem)
end

# `elsize`: part of the strided arrays interface, used for `pointer`

Base.elsize(::Type{A}) where {A<:FixedSizeArray} = Base.elsize(parent_type(A))

# `reshape`: specializing it to ensure it returns a `FixedSizeArray`

function Base.reshape(a::FixedSizeArray, size::(NTuple{N,Int} where {N}))
    len = checked_dims(size)
    if length(a) != len
        throw(DimensionMismatch("new shape not consistent with existing array length"))
    end
    new_fixed_size_array(a.mem, size)
end

# `iterate`: the `AbstractArray` fallback doesn't perform well, so add our own methods

function Base.iterate(a::FixedSizeArray)
    iterate(a.mem)
end
function Base.iterate(a::FixedSizeArray, state)
    iterate(a.mem, state)
end

const FixedSizeArrayDefault = FixedSizeArray{T, N, default_underlying_storage_type{T}} where {T, N}
const FixedSizeVectorDefault = FixedSizeArrayDefault{T, 1} where {T}
const FixedSizeMatrixDefault = FixedSizeArrayDefault{T, 2} where {T}
