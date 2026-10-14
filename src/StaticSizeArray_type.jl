# Dimensions of a `StaticSizeArray`: each one is either a `StaticInt`, a compile-time
# constant which takes no space, or an `Int`, e.g. `Tuple{StaticInt{4}, Int}` for `(4, n)`.
const MixedDims{N} = NTuple{N, Union{Int, StaticInt}}

"""
    StaticSizeArray{T,N,D,Mem<:DenseVector{T}} <: DenseArray{T,N}

An `N`-dimensional array with element type `T`, stored in column-major order in the memory
backend `Mem`, whose dimensions are given by a tuple `dims::D` of `Int`s and `StaticInt`s
from [Static.jl](https://github.com/SciML/Static.jl). Each dimension can be static, i.e. a
compile-time constant which is part of the type `D` and takes no space in the array, or
dynamic, i.e. an `Int` stored in the array. E.g. with `dims = (static(5), static(4), n)`,
`D = Tuple{StaticInt{5}, StaticInt{4}, Int}`, the array has size `(5, 4, n)`, and `n` can be
different for different arrays of the same type.

The static dimensions allow the compiler to compute the linear indices, the strides, and
the sizes in these dimensions at compile time, which is useful e.g. for arrays where some
dimensions are fixed (variables, nodes, ...) while others (elements of a mesh, particles,
...) change at runtime.

Construct a `StaticSizeArray` with [`FixedSizeArrays.wrap`](@ref) around an existing memory
backend (e.g. a `Memory` or a `Vector`), with an `undef` constructor, or from a
[`FixedSizeArray`](@ref):

    StaticSizeArray{T}(undef, dims::Tuple)
    StaticSizeArray{T,N,D,Mem}(undef, dims::Tuple)
    StaticSizeArray(a::FixedSizeArray, dims::Tuple)

The `undef` constructors use the default memory backend (`Memory{T}` on Julia v1.11+,
`Vector{T}` on Julia v1.10), unless `Mem` is given. `StaticSizeArray(a, dims)` returns a
`StaticSizeArray` which shares the memory backend with `a`, i.e. it does not copy the
elements. `size(a)` must be equal to `dims`, otherwise a `DimensionMismatch` is thrown.

!!! note

    `StaticSizeArray` is public, but not exported. Use it as `FixedSizeArrays.StaticSizeArray`,
    or import it with `using FixedSizeArrays: StaticSizeArray`.

# Examples

```jldoctest
julia> using FixedSizeArrays: FixedSizeArrays, FixedSizeArray, StaticSizeArray

julia> using Static: static

julia> a = FixedSizeArrays.wrap(StaticSizeArray, collect(1:12), (static(2), static(3), 2));

julia> size(a), a[2, 3, 2]
((2, 3, 2), 12)

julia> b = StaticSizeArray(FixedSizeArray([1.0 2.0; 3.0 4.0]), (2, static(2)));  # shares the memory

julia> b == [1.0 2.0; 3.0 4.0]
true
```
"""
# (`D` has no bound in the type, so that e.g. `StaticSizeArray{T, N, D, Mem} where {D}` is a
# `StaticSizeArray`; the inner constructor only accepts dimensions of type `MixedDims{N}`.)
struct StaticSizeArray{T,N,D,Mem<:DenseVector{T}} <: DenseArray{T,N}
    mem::Mem
    dims::D
    global function new_static_size_array(mem::DenseVector{T}, dims::MixedDims{N}) where {T,N}
        new{T,N,typeof(dims),typeof(mem)}(mem, dims)
    end
end
