# Usage

## The `FixedSizeArray` type

### Constructors

The basic type provided by this package is [`FixedSizeArray`](@ref).
Its constructor uses a similar syntax to `Base`'s `Array`:

```@repl
using FixedSizeArrays
v = FixedSizeArray{Float64}(undef, 3)
M = FixedSizeArray{Float64}(undef, 2, 2)
```

You can also construct `FixedSizeArray`s from other arrays:

```jldoctest
julia> using FixedSizeArrays

julia> arr = [10 20; 30 14]
2×2 Matrix{Int64}:
 10  20
 30  14

julia> FixedSizeArray(arr)  # construct from an `AbstractArray` value
2×2 FixedSizeArray{Int64, 2, Memory{Int64}}:
 10  20
 30  14

julia> FixedSizeArray{Float64}(arr)  # construct from an `AbstractArray` value while converting element type
2×2 FixedSizeArray{Float64, 2, Memory{Float64}}:
 10.0  20.0
 30.0  14.0
```

For 1- and 2-dimensional arrays you can use the aliases [`FixedSizeVector`](@ref) and [`FixedSizeMatrix`](@ref), respectively:

```@repl
using FixedSizeArrays
v = FixedSizeVector{Float64}(undef, 2)
M = FixedSizeMatrix{Float64}(undef, 3, 3)
```

### The memory backend

The `FixedSizeArray{T,N,Mem}` type has three parameters:

* `T`: the type of the elements of the array (e.g. `Int64` or `Float64`);
* `N`: the number of dimensions of the array (e.g. 1 for a vector, 2 for a matrix);
* `Mem<:DenseVector{T}`: the memory backend.

!!! note "Implementation details"

    In Julia v1.11+, the default memory backend is the [`Memory{T}`](https://docs.julialang.org/en/v1/base/arrays/#Core.Memory) type.
    Since Julia v1.10 does not have the `Memory` type, to make this package usable also on Julia v1.10 `Vector{T}` is used as memory backend, but many of the memory/performance optimizations enabled by this package will not be available in that version of Julia and in general `FixedSizeArrays.jl` does not provide significant improvements compared to `Base`'s `Array` for that specific version.

To make it easier to refer to the concrete type `FixedSizeArray{T,N,Mem}` with the default memory backend, the following convenient aliases are available:

* [`FixedSizeArrayDefault{T,N}`](@ref);
* [`FixedSizeVectorDefault{T}`](@ref);
* [`FixedSizeMatrixDefault{T}`](@ref).

## The `collect_as` utility

The [array literals syntax](https://docs.julialang.org/en/v1/manual/arrays/#man-array-literals) `[A, B, C, ...]` to construct arrays is limited to `Base`'s `Array` and cannot be extended to custom array types.
`FixedSizeArrays.jl` provides a convenient function [`collect_as`](@ref) to overcome this limitation and construct `FixedSizeArray`s out of any iterable:

```jldoctest
julia> iter = (i for i ∈ 7:9 if i≠8);

julia> using FixedSizeArrays

julia> const ca = FixedSizeArrays.collect_as;

julia> ca(FixedSizeArray, iter)  # construct from an arbitrary iterator
2-element FixedSizeArray{Int64, 1, Memory{Int64}}:
 7
 9

julia> ca(FixedSizeArray{Float64}, iter)  # construct from an arbitrary iterator while converting element type
2-element FixedSizeArray{Float64, 1, Memory{Float64}}:
 7.0
 9.0

julia> ca(FixedSizeVectorDefault, (3.14, -4.2, 2.68))  # construct from a tuple
3-element FixedSizeArray{Float64, 1, Memory{Float64}}:
  3.14
 -4.2
  2.68
```

## `BoundsErrorLight` exception

To facilitate the [escape analysis](https://en.wikipedia.org/wiki/Escape_analysis) of `FixedSizeArray`s, accessing an out-of-bound index of these arrays raises a [`BoundsErrorLight`](@ref) exception when possible.
This exception type does not store the entire array for reporting the error message, thus enabling more performance optimizations compared to arrays which throw [`BoundsError`](https://docs.julialang.org/en/v1/base/base/#Core.BoundsError) exceptions.

!!! warning "Compatibility of exception raised"

    We do not guarantee to throw either `BoundsError` or `BoundsErrorLight` when accessing an out-of-bound index of a `FixedSizeArray`.
    The exact exception raised may change at any point without breaking the semantic versioning compatibility contract.
