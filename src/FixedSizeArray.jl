export
    FixedSizeArray, FixedSizeVector, FixedSizeMatrix,
    FixedSizeArrayDefault, FixedSizeVectorDefault, FixedSizeMatrixDefault

struct FixedSizeArray{T,N,Mem<:DenseVector{T}} <: DenseArray{T,N}
    mem::Mem
    size::NTuple{N,Int}
    global function new_fixed_size_array(mem::DenseVector{T}, size::NTuple{N,Int}) where {T,N}
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
            Type{FixedSizeArray{T, N, Mem{T}} where {N}} where {T},
            Type{FixedSizeArray{T, N, Mem{T}} where {T}} where {N},
            # 2 fixed parameters
            Type{FixedSizeArray{T, N, Mem{T}}} where {T, N},
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

function check_constructor_is_allowed(::Type{T}) where {T <: FixedSizeArray}
    if Type{T} <: FixedSizeArrayAllowedConstructorType
        T
    else
        throw(ArgumentError("technical limitation: type is not among the allowed `FixedSizeArray` constructors (try deleting type parameter constraints)"))
    end
end

const FixedSizeVector{T} = FixedSizeArray{T,1}
const FixedSizeMatrix{T} = FixedSizeArray{T,2}

function parent_type_with_default(::Type{<:(FixedSizeArray{E, N, T} where {N})}) where {E, T <: DenseVector{E}}
    T
end
function parent_type_with_default(::Type{<:FixedSizeArray{E}}) where {E}
    default_underlying_storage_type{E}
end
function parent_type_with_default(::Type{<:FixedSizeArray})
    default_underlying_storage_type
end
for T ∈ (Vector, optional_memory...)
    FSA = FixedSizeArray{E, N, T{E}} where {E, N}
    t_fsa = Union{
        Type{FSA},
        Type{FSA{E, N} where {E}} where {N},
    }
    @eval begin
        function parent_type_with_default(::$t_fsa)
            $T
        end
    end
end

function check_ndims(::Type{<:(FixedSizeArray{E, N} where {E})}, size::Tuple{Vararg{Integer}}) where {N}
    if size isa Tuple{Vararg{Any, N}}
        size
    else
        throw(DimensionMismatch("mismatch between dimension count in type and the length of the size tuple"))
    end
end
function check_ndims(::Type{<:FixedSizeArray}, size::Tuple{Vararg{Integer}})
    size
end

function undef_constructor(::Type{FSA}, size::Tuple{Vararg{Integer}}) where {E, FSA <: FixedSizeArray{E}}
    size = check_ndims(FSA, size)
    s = map(Int, size)
    Mem = parent_type_with_default(FSA)
    len = checked_dims(s)
    mem = Mem(undef, len)
    new_fixed_size_array(mem, s)
end

function (::Type{FSA})(::UndefInitializer, size::Tuple{Vararg{Integer}}) where {E, FSA <: FixedSizeArray{E}}
    undef_constructor(FSA, size)
end
function (::Type{FSA})(::UndefInitializer, size::Vararg{Integer}) where {E, FSA <: FixedSizeArray{E}}
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

function (::Type{FSA})(src::AbstractArray) where {FSA <: FixedSizeArray}
    axs = axes(src)
    if !axes_are_one_based(axs)
        throw(DimensionMismatch("source array has a non-one-based indexing axis"))
    end
    collect_as(FSA, src)
end

# conversion

Base.convert(::Type{T}, a::T) where {T<:FixedSizeArray} = a
Base.convert(::Type{T}, a::AbstractArray) where {T<:FixedSizeArray} = T(a)

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
