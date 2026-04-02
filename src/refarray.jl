# Internal type used as the default, heap-backed parent container for FixedSizeArray.
# This is used over Memory{T}, because a MemoryRef inlines into the parent struct, and
# contains the pointer. This way, it saves a memory indirection over Memory.
# This type does not boundscheck, hence the unsafety. The FSA wrapper will handle boundschecks.
struct UnsafeRefArray{T} <: DenseVector{T}
    ref::MemoryRef{T}
end

# TODO: Would be nice to be able to omit the check for n < 0 in the Memory
# constructor, but this is in the memorynew built-in and cannot be disabled
function UnsafeRefArray{T}(::UndefInitializer, n::Int) where T
    mem = Memory{T}(undef, n)
    UnsafeRefArray{T}(memoryref(mem))
end

UnsafeRefArray{T}(::UndefInitializer, n::Tuple{Int}) where T = UnsafeRefArray{T}(undef, @inbounds(n[1]))

@static if VERSION < v"1.12.0-DEV.966"
    Base.parent(x::UnsafeRefArray) = x.ref.mem
else
    Base.parent(x::UnsafeRefArray) = parent(x.ref)
end

# This version introduced Base.memoryindex as public API
@static if VERSION < v"1.13.0-DEV.1289"
    internal_memindex(x::MemoryRef) = Core.memoryrefoffset(x)
else
    internal_memindex(x::MemoryRef) = Base.memoryindex(x)
end

function Base.checkbounds(A::UnsafeRefArray, is...)
    checkbounds_lightboundserror(A, is...)
end

Base.size(x::UnsafeRefArray) = ((length(parent(x)) - internal_memindex(x.ref) + 1),)

Base.firstindex(::UnsafeRefArray) = 1

function Base.getindex(x::UnsafeRefArray, i::Int)
    @inbounds(memoryref(x.ref, i)[])
end

function Base.setindex!(x::UnsafeRefArray{T}, v, i::Int) where T
    vT = convert(T, v)::T
    @inbounds(memoryref(x.ref, i)[] = vT)
end

Base.IndexStyle(::Type{<:UnsafeRefArray}) = Base.IndexLinear()

Base.elsize(::Type{UnsafeRefArray{T}}) where {T} = Base.elsize(Memory{T})

function Base.similar(A::UnsafeRefArray, ::Type{S}, dims::Dims) where S
    UnsafeRefArray{S}(undef, dims)
end

function Base.isassigned(x::UnsafeRefArray, i::Int)
    @boundscheck checkbounds(Bool, x, i) || false
    ref = @inbounds memoryref(x.ref, i)
    isassigned(ref)
end

function Base.copy(x::UnsafeRefArray)
    len = length(x)
    newmem = Memory{eltype(x)}(undef, len)
    newref = memoryref(newmem)
    @inbounds unsafe_copyto!(newref, x.ref, len)
    typeof(x)(newref)
end

function Base.unsafe_convert(::Type{Ptr{T}}, x::UnsafeRefArray) where T
    Base.unsafe_convert(Ptr{T}, x.ref)
end

Base.dataids(x::UnsafeRefArray) = Base.dataids(parent(x))

# Collects.jl interface for UnsafeRefArray

function (c::Collect)(::Type{UnsafeRefArray{T}}, collection) where {T}
    if collection isa UnsafeRefArray{T}
        return copy(collection)
    end
    if (T isa Type) && Base.IteratorSize(collection) isa Union{Base.HasLength, Base.HasShape}
        result = UnsafeRefArray{T}(undef, Int(length(collection))::Int)
        copyto!(result, collection)
        return result
    end
    vec = c(Vector{T}, collection)
    result = UnsafeRefArray{T}(undef, length(vec))
    copyto!(result, vec)
    result
end

function (c::Collect)(::Type{UnsafeRefArray}, collection)
    if collection isa UnsafeRefArray
        return copy(collection)
    end
    vec = c(Vector, collection)
    T = eltype(vec)
    result = UnsafeRefArray{T}(undef, length(vec))
    copyto!(result, vec)
    result
end
