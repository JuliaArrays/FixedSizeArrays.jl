module FixedSizeArraysGPUArraysCoreExt

# Support for `FixedSizeArray`s whose memory backend is a GPU array, i.e. a
# `GPUArraysCore.AbstractGPUArray` (e.g. a `CuVector` or a `JLVector`).  The generic
# `AbstractArray` fallbacks of many operations (broadcasting, `fill!`, reductions, `show`,
# ...) index the array element by element, which is not allowed (or extremely slow) for
# GPU arrays.  Instead, forward these operations to the parent GPU vector, or to the
# `N`-dimensional GPU array obtained by reshaping it, for which the GPU array packages
# implement them with kernels.

using Adapt: Adapt, adapt
using GPUArraysCore: AbstractGPUArray, AbstractGPUArrayStyle
using FixedSizeArrays: FixedSizeArrays, FixedSizeArray, checked_dims, wrap
using Base.Broadcast: Broadcast, Broadcasted, BroadcastStyle

const GPUFixedSizeArray{T, N} = FixedSizeArray{T, N, <:AbstractGPUArray}

# The `N`-dimensional GPU array sharing the memory of `a`.
gpu_array(a::GPUFixedSizeArray) = reshape(parent(a), size(a))

# Wrap the result of an operation on GPU arrays again, if it is an array.
rewrap(r::AbstractGPUArray) = wrap(FixedSizeArray, vec(r), size(r))
rewrap(r) = r

# Adapt: adapt the memory backend and rewrap it, so that e.g. converting a GPU-backed
# `FixedSizeArray` for use in a kernel (`CUDA.cudaconvert`, KernelAbstractions.jl, ...)
# gives a `FixedSizeArray` backed by a device array, and `adapt(Array, a)` gives a
# `FixedSizeArray` backed by a `Vector`.  This is restricted to GPU-backed arrays: for
# other memory backends, `Adapt.adapt` keeps treating a `FixedSizeArray` as an array to
# convert as a whole, e.g. `adapt(CuArray, a)` gives a `CuArray`.

function Adapt.adapt_structure(to, a::GPUFixedSizeArray)
    wrap(FixedSizeArray, adapt(to, parent(a)), size(a))
end

# broadcasting: a GPU broadcast style wrapping the style of the memory backend, so that
# broadcasting runs as a kernel of the memory backend, while out-of-place broadcasting
# returns a `FixedSizeArray`

struct FixedSizeArrayGPUBroadcastStyle{N, S <: AbstractGPUArrayStyle{N}} <: AbstractGPUArrayStyle{N} end

function (::Type{<:(FixedSizeArrayGPUBroadcastStyle{N, S} where {N})})(::Val{M}) where {S, M}
    s = S(Val(M))
    FixedSizeArrayGPUBroadcastStyle{M, typeof(s)}()
end

wrap_style(s::AbstractGPUArrayStyle{N}) where {N} = FixedSizeArrayGPUBroadcastStyle{N, typeof(s)}()
wrap_style(s::FixedSizeArrayGPUBroadcastStyle) = s
wrap_style(s) = s

memory_style(::FixedSizeArrayGPUBroadcastStyle{N, S}) where {N, S} = S()

function Base.BroadcastStyle(::Union{Type{<:FixedSizeArray{T, N, Mem}},
                                     Type{<:SubArray{T, N, FixedSizeArray{T, N, Mem}}}}) where {T, N, Mem <: AbstractGPUArray}
    s = BroadcastStyle(Mem)
    wrap_style(typeof(s)(Val(N)))
end

# Combine with other GPU styles via the rules of the memory backend.
function Base.BroadcastStyle(a::FixedSizeArrayGPUBroadcastStyle, b::AbstractGPUArrayStyle)
    wrap_style(Broadcast.result_style(memory_style(a), b))
end
function Base.BroadcastStyle(a::FixedSizeArrayGPUBroadcastStyle, b::FixedSizeArrayGPUBroadcastStyle)
    wrap_style(Broadcast.result_style(memory_style(a), memory_style(b)))
end
# resolve the ambiguity with `BroadcastStyle(::A, ::A) where {A <: AbstractArrayStyle}`
Base.BroadcastStyle(a::S, ::S) where {S <: FixedSizeArrayGPUBroadcastStyle} = a

function Base.similar(bc::Broadcasted{<:FixedSizeArrayGPUBroadcastStyle}, ::Type{T}) where {T}
    similar(bc, T, axes(bc))
end
function Base.similar(bc::Broadcasted{FixedSizeArrayGPUBroadcastStyle{N, S}}, ::Type{T}, dims) where {N, S, T}
    bc_memory = Broadcasted{S}(bc.f, bc.args, bc.axes)
    rewrap(similar(bc_memory, T, dims))
end

# `similar`: the generic method for `FixedSizeArray`s would construct the memory backend
# with an `undef` constructor, which GPU array types may not have.

function Base.similar(a::GPUFixedSizeArray, ::Type{E}, size::NTuple{N, Int}) where {E, N}
    wrap(FixedSizeArray, similar(parent(a), E, (checked_dims(size),)), size)
end

# `fill!` and `copyto!`

function Base.fill!(a::GPUFixedSizeArray, x)
    fill!(parent(a), x)
    a
end

function Base.copyto!(dst::GPUFixedSizeArray, src::AbstractGPUArray)
    copyto!(parent(dst), src)
    dst
end
function Base.copyto!(dst::AbstractGPUArray, src::GPUFixedSizeArray)
    copyto!(dst, parent(src))
end
function Base.copyto!(dst::GPUFixedSizeArray, doff::Integer, src::AbstractGPUArray, soff::Integer, n::Integer)
    copyto!(parent(dst), doff, src, soff, n)
    dst
end
function Base.copyto!(dst::AbstractGPUArray, doff::Integer, src::GPUFixedSizeArray, soff::Integer, n::Integer)
    copyto!(dst, doff, parent(src), soff, n)
end

# `map`, reductions and comparisons

Base.map(f, a::GPUFixedSizeArray) = Broadcast.broadcast_preserving_zero_d(f, a)

function Base.mapreduce(f, op, a::GPUFixedSizeArray; kwargs...)
    rewrap(mapreduce(f, op, gpu_array(a); kwargs...))
end

Base.any(a::GPUFixedSizeArray{Bool}) = any(gpu_array(a))
Base.all(a::GPUFixedSizeArray{Bool}) = all(gpu_array(a))
Base.any(f::Function, a::GPUFixedSizeArray) = any(f, gpu_array(a))
Base.all(f::Function, a::GPUFixedSizeArray) = all(f, gpu_array(a))
function Base.count(f::Function, a::GPUFixedSizeArray; kwargs...)
    rewrap(count(f, gpu_array(a); kwargs...))
end

gpu_array(a::AbstractGPUArray) = a

for (A, B) ∈ ((GPUFixedSizeArray, GPUFixedSizeArray), (GPUFixedSizeArray, AbstractGPUArray), (AbstractGPUArray, GPUFixedSizeArray))
    @eval begin
        function Base.:(==)(a::$A, b::$B)
            (axes(a) == axes(b)) && (gpu_array(a) == gpu_array(b))
        end
        function Base.isapprox(a::$A, b::$B; kwargs...)
            (axes(a) == axes(b)) && isapprox(gpu_array(a), gpu_array(b); kwargs...)
        end
    end
end

# `collect` and `show`: copy to the host first (`Array(a)` already works, via `copyto!`)

Base.collect(a::GPUFixedSizeArray{T, N}) where {T, N} = copyto!(Array{T, N}(undef, size(a)), a)

host_array(a::GPUFixedSizeArray) = wrap(FixedSizeArray, Vector(parent(a)), size(a))

Base.print_array(io::IO, a::GPUFixedSizeArray) = Base.print_array(io, host_array(a))
Base._show_nonempty(io::IO, a::GPUFixedSizeArray, prefix::String) = Base._show_nonempty(io, host_array(a), prefix)
Base._show_empty(io::IO, a::GPUFixedSizeArray) = Base._show_empty(io, host_array(a))
Base.show_vector(io::IO, a::GPUFixedSizeArray, args...) = Base.show_vector(io, host_array(a), args...)

end # module FixedSizeArraysGPUArraysCoreExt
