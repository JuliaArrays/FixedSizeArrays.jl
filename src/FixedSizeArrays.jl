module FixedSizeArrays

export FixedSizeArray, FixedSizeVector, FixedSizeMatrix

struct FixedSizeArray{T,N} <: DenseArray{T,N}
    mem::Memory{T}
    size::NTuple{N,Int}
    function FixedSizeArray{T,N}(mem::Memory{T}, size::NTuple{N,Int}) where {T,N}
        new{T,N}(mem, size)
    end
end

const FixedSizeVector{T} = FixedSizeArray{T,1}
const FixedSizeMatrix{T} = FixedSizeArray{T,2}

function FixedSizeArray{T,N}(::UndefInitializer, size::NTuple{N,Int}) where {T,N}
    FixedSizeArray{T,N}(Memory{T}(undef, prod(size)), size)
end
function FixedSizeArray{T,N}(::UndefInitializer, size::Vararg{Int,N}) where {T,N}
    FixedSizeArray{T,N}(undef, size)
end
function FixedSizeArray{T}(::UndefInitializer, size::Vararg{Int,N}) where {T,N}
    FixedSizeArray{T,N}(undef, size)
end
function FixedSizeArray{T}(::UndefInitializer, size::NTuple{N,Int}) where {T,N}
    FixedSizeArray{T,N}(undef, size)
end

Base.IndexStyle(::Type{<:FixedSizeArray}) = IndexLinear()
Base.@propagate_inbounds Base.getindex(A::FixedSizeArray, i::Int) = A.mem[i]
Base.@propagate_inbounds Base.setindex!(A::FixedSizeArray, v, i::Int) = A.mem[i] = v

Base.size(a::FixedSizeArray) = getfield(a, :size)

function Base.similar(::FixedSizeArray, ::Type{S}, size::NTuple{N,Int}) where {S,N}
    FixedSizeArray{S,N}(undef, size...)
end

end # module FixedSizeArrays
