# `GenericMemoryRef` memory backend (Julia v1.11 and later)
#
# A `FixedSizeArrayDims{T,N,D,<:GenericMemoryRef{kind,T}}` (e.g. a
# `FixedSizeArray{T,N,MemoryRef{T}}`) stores a reference to its first element: the `Memory`
# holding the elements, and a pointer (or offset) into it. Element `i` of the array is
# element `i` after the reference. The dimensions of the array give the length, there is
# no extra length field.
#
# Invariant (established by all constructors, relied upon by indexing without bounds
# checks): `length(a)` elements starting at the reference are within the memory, i.e.
# `memoryrefoffset(ref) - 1 + length(a) <= length(ref.mem)`.
#
# Unlike a `Vector` backend, which is a mutable object whose `ref` field `resize!` can
# change, the reference is an immutable field of the immutable array: the data
# pointer is loop invariant, and the `Memory` it points into stays alive (rooted) while
# the array is.

const memoryrefnew = @static if isdefined(Core, :memoryrefnew)
    Core.memoryrefnew
else
    Core.memoryref
end

const MemRef{T} = GenericMemoryRef{K, T, AS} where {K, AS}

# The `GenericMemory` type corresponding to a `GenericMemoryRef` type. The element type
# may be a type variable.
memory_type(::Type{GenericMemoryRef{K, T, AS}}) where {K, T, AS} = GenericMemory{K, T, AS}
memory_type(::Type{<:(GenericMemoryRef{K, T, AS} where {T})}) where {K, AS} = (GenericMemory{K, T, AS} where {T})

# Number of elements of the memory from the reference (inclusive) to its end.
function available_length(ref::GenericMemoryRef)
    length(getfield(ref, :mem)) - Base.memoryrefoffset(ref) + 1
end

@inline function element_ref(ref::GenericMemoryRef, i::Int)
    memoryrefnew(ref, i, false)
end

@inline function storage_getindex(ref::GenericMemoryRef, i::Int)
    Core.memoryrefget(element_ref(ref, i), Base.default_access_order(ref), false)
end
@inline function storage_setindex!(ref::MemRef{T}, x, i::Int) where {T}
    val = x isa T ? x : convert(T, x)::T
    Core.memoryrefset!(element_ref(ref, i), val, Base.default_access_order(ref), false)
    nothing
end
@inline function storage_isassigned(ref::GenericMemoryRef, i::Int)
    Core.memoryref_isassigned(element_ref(ref, i), Base.default_access_order(ref), false)
end

function new_storage(::Type{Mem}, len::Int) where {Mem <: GenericMemoryRef}
    memoryrefnew(memory_type(Mem)(undef, len))
end

function similar_storage(::GenericMemoryRef{K, <:Any, AS}, ::Type{E}, len::Int) where {K, AS, E}
    memoryrefnew(GenericMemory{K, E, AS}(undef, len))
end

# A fresh memory of exactly `len` elements, starting at offset 1.
function copy_storage(ref::GenericMemoryRef, len::Int)
    mem = getfield(ref, :mem)
    new = memoryrefnew(typeof(mem)(undef, len))
    if iszero(len)
    elseif new isa MemoryRef
        unsafe_copyto!(new, ref, len)
    else  # e.g. atomic memory: copy element by element
        for i ∈ 1:len
            if storage_isassigned(ref, i)
                storage_setindex!(new, storage_getindex(ref, i), i)
            end
        end
    end
    new
end

function collect_storage(collect, ::Type{Mem}, iterator) where {Mem <: GenericMemoryRef}
    memoryrefnew(collect(memory_type(Mem), iterator))
end

storage_elsize(::Type{<:GenericMemoryRef{K, T, AS}}) where {K, T, AS} = Base.elsize(GenericMemory{K, T, AS})

storage_array_type(::Type{R}) where {R <: GenericMemoryRef} = memory_type(R)


storage_dataids(ref::GenericMemoryRef) = Base.dataids(getfield(ref, :mem))

# `parent(a)` is the reference, so the indices into the parent are `1:length(a)`
# (relative to the reference).
storage_parentindices(::GenericMemoryRef, f::FixedSizeArrayDims) = (Base.OneTo(length(f)),)

# A view of the elements of the memory, for generic fallbacks (`copyto!`, `rand!`, ...).
function storage_vector(ref::GenericMemoryRef, len::Int)
    off = Base.memoryrefoffset(ref)
    view(getfield(ref, :mem), off:(off + len - 1))
end

function with_stripped_type_parameters_unchecked(::TypeParametersElementType, ::Type{<:(GenericMemoryRef{K, T, AS} where {T})}) where {K, AS}
    s = GenericMemoryRef{K, T, AS} where {T}
    Val{s}()
end

# `copyto!` between memory references: `memmove` for `isbits` element types, otherwise
# with the necessary write barriers, overlapping ranges are handled.
const MemRefFixedSizeArray{T} = FixedSizeArrayDims{T, N, D, <:MemRef{T}} where {N, D}

function maybe_memoryref(a::FixedSizeArrayDims)
    p = parent(a)
    if p isa GenericMemoryRef
        p
    else
        maybe_memoryref(p)
    end
end
maybe_memoryref(a::Array) = getfield(a, :ref)
maybe_memoryref(m::GenericMemory) = memoryrefnew(m)
maybe_memoryref(::Any) = nothing

@noinline function throw_copyto_negative(n)
    throw(ArgumentError(lazy"tried to copy n=$n elements, but n should be non-negative"))
end

Base.@propagate_inbounds function copyto_memoryref!(dst, doff::Int, src, soff::Int, n::Int)
    n < 0 && throw_copyto_negative(n)
    iszero(n) && return dst
    @boundscheck checkbounds(dst, doff:(doff + n - 1))
    @boundscheck checkbounds(src, soff:(soff + n - 1))
    dr = maybe_memoryref(dst)
    sr = maybe_memoryref(src)
    unsafe_copyto!(element_ref(dr, doff), element_ref(sr, soff), n)
    dst
end

# use the fast path only if one side has a `GenericMemoryRef` backend and both sides are
# backed by a memory reference with the same element type and memory kind
function copyto_fast_path(dst, src)
    dr = maybe_memoryref(dst)
    sr = maybe_memoryref(src)
    (dr isa MemoryRef) && (sr isa MemoryRef) && (eltype(dr) === eltype(sr))
end

use_memoryref_copy(dst::MemRefFixedSizeArray, src) = copyto_fast_path(dst, src)
use_memoryref_copy(dst, src::MemRefFixedSizeArray) = copyto_fast_path(dst, src)
use_memoryref_copy(dst::MemRefFixedSizeArray, src::MemRefFixedSizeArray) = copyto_fast_path(dst, src)

# `wrap`

@noinline function throw_wrap_ref_too_short(avail, s)
    throw(DimensionMismatch(lazy"the memory has $avail elements from the reference on, fewer than the product of the dimensions $s"))
end

# A `Colon` dimension is computed from all elements of the memory from `ref` on.
function wrap(::Type{FSA}, ref::GenericMemoryRef, dims::Tuple) where {FSA <: FixedSizeArrayDims}
    avail = available_length(ref)
    d = convert_dims(dims_type(FSA), resolve_colon(check_ndims(FSA, dims), avail))
    s = map(dim_value, d)
    len = checked_dims(s)
    if avail < len
        throw_wrap_ref_too_short(avail, s)
    end
    check_array_type(FSA, new_fixed_size_array(ref, d))
end
function wrap(::Type{FSA}, ref::GenericMemoryRef, size1::Union{Integer, Val}, size::Vararg{Union{Integer, Val}}) where {FSA <: FixedSizeArrayDims}
    wrap(FSA, ref, (size1, size...))
end

# A `Vector` is wrapped with a `GenericMemoryRef` backend (its current reference) when the
# requested type has such a backend, e.g. `FixedSizeArray{T, N, MemoryRef{T}} where {T, N}`
# or `FixedSizeArrayDims{T, N, D, MemoryRef{T}} where {T, N, D}`.
function wrap_backend(::Type{<:(FixedSizeArrayDims{T, N, D, <:MemRef{T}} where {T, N, D})}, parent::Vector)
    getfield(parent, :ref)
end
