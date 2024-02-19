module FixedSizeArrays

export FixedSizeArray, FixedSizeVector, FixedSizeMatrix

mutable struct FixedSizeArray{T,N} <: DenseArray{T,N}
    ref::MemoryRef{T}
    const size::NTuple{N,Int}
end

const FixedSizeVector{T} = FixedSizeArray{T,1}
const FixedSizeMatrix{T} = FixedSizeArray{T,2}

eval(:(function (self::Type{FixedSizeArray{T,N}})(::UndefInitializer, size::Vararg{Int,N}) where {T,N}
    mem = fieldtype(fieldtype(self, :ref), :mem)(undef, prod(size))
    return $(Expr(:new, :self, :(Core.memoryref(mem)), :(size)))
end))

function Base.setindex!(A::FixedSizeArray{T}, x, i::Int) where {T}
    Base.@_noub_if_noinbounds_meta
    @boundscheck (i - 1)%UInt < length(A)%UInt || throw_boundserror(A, (i,))
    Core.memoryrefset!(Core.memoryref(A.ref, i, false), x isa T ? x : convert(T,x)::T, :not_atomic, false)
    return A
end
function Base.setindex!(A::FixedSizeArray{T}, x, i1::Int, i2::Int, I::Int...) where {T}
    @inline
    Base.@_noub_if_noinbounds_meta
    @boundscheck checkbounds(A, i1, i2, I...) # generally _to_linear_index requires bounds checking
    Core.memoryrefset!(Core.memoryref(A.ref, Base._to_linear_index(A, i1, i2, I...), false), x isa T ? x : convert(T,x)::T, :not_atomic, false)
    return A
end

function Base.getindex(A::FixedSizeArray, i::Int)
    Base.@_noub_if_noinbounds_meta
    @boundscheck Base.ult_int(Base.bitcast(UInt, Base.sub_int(i, 1)), Base.bitcast(UInt, length(A))) || throw_boundserror(A, (i,))
    Core.memoryrefget(Core.memoryref(getfield(A, :ref), i, false), :not_atomic, false)
end
function Base.getindex(A::FixedSizeArray, i1::Int, i2::Int, I::Int...)
    @inline
    @boundscheck checkbounds(A, i1, i2, I...) # generally _to_linear_index requires bounds checking
    return @inbounds A[Base._to_linear_index(A, i1, i2, I...)]
end

Base.size(a::FixedSizeArray) = getfield(a, :size)

end # module FixedSizeArrays
