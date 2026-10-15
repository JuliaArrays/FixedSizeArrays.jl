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

To "convert" a wrapper type like `Adjoint` while preserving the wrapper structure, you can use `Adapt`:

```jldoctest
julia> using FixedSizeArrays, Adapt

julia> arr = transpose([10 20; 30 14])
2×2 transpose(::Matrix{Int64}) with eltype Int64:
 10  30
 20  14

julia> adapt(FixedSizeArray, arr)
2×2 transpose(::FixedSizeArray{Int64, 2, Memory{Int64}}) with eltype Int64:
 10  30
 20  14
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

## Static dimensions: the `FixedSizeArrayDims` type

Many arrays in numerical codes have a shape that is partly known in advance.
Take the solution of a discontinuous Galerkin method in three dimensions, stored as `u[variable, i, j, k, element]`: the number of variables (five for the compressible Euler equations) and of nodes per direction (say four) are fixed by the choice of the equations and the polynomial degree, but the number of elements depends on the mesh, and may change while the simulation runs, e.g. with adaptive mesh refinement.
With a plain `Array` or `FixedSizeArray` the compiler knows none of these sizes, so every index computation multiplies with sizes loaded from memory, and loops over the small dimensions cannot be fully unrolled.

For such arrays, the dimensions of the underlying type [`FixedSizeArrays.FixedSizeArrayDims{T,N,D,Mem}`](@ref FixedSizeArrays.FixedSizeArrayDims) can be a tuple of type `D` of `Int`s and `StaticInt`s from [Static.jl](https://github.com/SciML/Static.jl).
Any dimension can be static or dynamic.
`FixedSizeArray{T,N,Mem}` is the alias `FixedSizeArrayDims{T,N,NTuple{N,Int},Mem}` with all dimensions dynamic, so all `FixedSizeArray`s are `FixedSizeArrayDims`s, but the arrays with static dimensions are not `FixedSizeArray`s.
`FixedSizeArrayDims` is public, but not exported.
Static.jl is a weak dependency of FixedSizeArrays.jl, so it is only loaded, and static dimensions are only available, once you load it yourself, e.g. with `using Static: static`.
The static dimensions are part of the type, and take no space in the array, so the compiler can compute the sizes in these dimensions, and the strides and linear indices involving them, at compile time:

```jldoctest
julia> using FixedSizeArrays: FixedSizeArrays, FixedSizeArrayDims

julia> using Static: static, StaticInt

julia> nelements = 3;

julia> u = FixedSizeArrayDims{Float64}(undef, (static(5), static(4), static(4), static(4), nelements));

julia> size(u), strides(u)
((5, 4, 4, 4, 3), (1, 5, 20, 80, 320))

julia> size(u, 2)  # a compile-time constant
4

julia> typeof(u)
FixedSizeArrayDims{Float64, 5, Tuple{StaticInt{5}, StaticInt{4}, StaticInt{4}, StaticInt{4}, Int64}, Memory{Float64}}
```

An array with static dimensions can also wrap an existing memory backend with [`FixedSizeArrays.wrap`](@ref), without copying.
At most one dimension can be given as `:`, which is then computed from the length of the memory backend, e.g. to wrap a `Vector` again after it was resized:

```jldoctest
julia> using FixedSizeArrays: FixedSizeArrays, FixedSizeArrayDims

julia> using Static: static

julia> v = collect(1.0:12.0);

julia> a = FixedSizeArrays.wrap(FixedSizeArrayDims, v, (static(2), static(3), :));

julia> size(a), parent(a) === v
((2, 3, 2), true)

julia> b = FixedSizeArrayDims{Int}(undef, (5, static(3)));  # a static dimension after a dynamic one

julia> size(b)
(5, 3)
```

The dimensions given to a constructor are converted to the dimensions type of the array, if its type gives it: e.g. `FixedSizeArray{Float64}(undef, (static(2), 3))` has two dynamic dimensions, while `FixedSizeArrayDims{Float64}(undef, (static(2), 3))` has a static one.
`similar(a)`, `similar(a, T)` and `similar(typeof(a), size)` keep the static dimensions, while `similar(a, T, size)` returns a `FixedSizeArray` with all dimensions dynamic, so that its type does not depend on the value of `size`.
Out-of-place broadcasting, `reshape` and `vec` return `FixedSizeArray`s, too.
With the package extension for [StaticArrayInterface.jl](https://github.com/JuliaArrays/StaticArrayInterface.jl), the static sizes and strides are available to packages such as LoopVectorization.jl, StrideArraysCore.jl and Polyester.jl.
Unlike for `FixedSizeArray`s, `Adapt.adapt` adapts the memory backend of an array with static dimensions and wraps it again, e.g. `adapt(CuArray, a)` returns a `FixedSizeArrayDims` backed by a `CuVector`.

### Comparison with other array types

Several packages in the Julia ecosystem encode sizes in the type.
They differ in which dimensions can be static, and in where the elements live:

| Array type | Static dimensions | Storage |
|:-----------|:------------------|:--------|
| `Array` | none, and a `Vector` can be resized | heap |
| [`FixedSizeArray`](@ref) | none, but the size never changes | any `DenseVector` memory backend |
| `SArray`, `MArray` ([StaticArrays.jl](https://github.com/JuliaArrays/StaticArrays.jl)) | all | inline (`SArray`, immutable) or a mutable object (`MArray`) |
| `HybridArray` ([HybridArrays.jl](https://github.com/JuliaArrays/HybridArrays.jl)) | any, using `StaticArrays.Dynamic()` for the dynamic ones | any parent `AbstractArray` |
| `StrideArray`, `PtrArray` ([StrideArraysCore.jl](https://github.com/JuliaSIMD/StrideArraysCore.jl)) | any, using `StaticInt`s; also static strides and offsets | a pointer (`PtrArray`, which must be kept alive with `GC.@preserve`), plus its memory for `StrideArray` |
| `FixedSizeArrayDims` | any, using `StaticInt`s | any `DenseVector` memory backend |

* StaticArrays.jl targets small arrays whose whole size is known, e.g. points or small matrices: the size is part of the type, and operations are unrolled.
  For large arrays, or ones with a dimension only known at runtime, it is not suitable.
* HybridArrays.jl also mixes static and dynamic dimensions, building on StaticArrays.jl, so that e.g. indexing along the static dimensions can return static arrays.
  It wraps any parent array, so it can be non-contiguous, while a `FixedSizeArrayDims` is always a `DenseArray` backed by the whole memory backend.
* StrideArraysCore.jl uses the same `Int`/`StaticInt` tuples for the sizes, and adds static strides and offsets, for LoopVectorization.jl.
  Its arrays are built around a pointer, while a `FixedSizeArrayDims` holds a reference to its memory backend, so the garbage collector knows about it, and the memory backend can be a GPU vector.
  The StaticArrayInterface.jl package extension makes `FixedSizeArrayDims`s work with StrideArraysCore.jl and LoopVectorization.jl, too.
* Compared to a `FixedSizeArray`, an array with static dimensions trades generality for speed: the static dimensions must be known at compile time, and each combination of them is a different type, so the code using it is compiled for each one.

## The `collect_as` utility

The [array literals syntax](https://docs.julialang.org/en/v1/manual/arrays/#man-array-literals) `[A, B, C, ...]` to construct arrays is limited to `Base`'s `Array` and cannot be extended to custom array types.
The package [Collects.jl](https://github.com/JuliaCollections/Collects.jl) provides a convenient function, `collect_as`, to overcome this limitation and construct `FixedSizeArray`s out of any iterable:

```jldoctest
julia> iter = (i for i ∈ 7:9 if i≠8);

julia> using FixedSizeArrays, Collects

julia> collect_as(FixedSizeArray, iter)  # construct from an arbitrary iterator
2-element FixedSizeArray{Int64, 1, Memory{Int64}}:
 7
 9

julia> collect_as(FixedSizeArray{Float64}, iter)  # construct from an arbitrary iterator while converting element type
2-element FixedSizeArray{Float64, 1, Memory{Float64}}:
 7.0
 9.0

julia> collect_as(FixedSizeVectorDefault, (3.14, -4.2, 2.68))  # construct from a tuple
3-element FixedSizeArray{Float64, 1, Memory{Float64}}:
  3.14
 -4.2
  2.68
```

See the Collects.jl Readme for more information.

## `BoundsErrorLight` exception

To facilitate the [escape analysis](https://en.wikipedia.org/wiki/Escape_analysis) of `FixedSizeArray`s, accessing an out-of-bound index of these arrays raises a `BoundsErrorLight` exception when possible.
This exception type does not store the entire array for reporting the error message, thus enabling more performance optimizations compared to arrays which throw [`BoundsError`](https://docs.julialang.org/en/v1/base/base/#Core.BoundsError) exceptions.

!!! warning "Compatibility of exception raised"

    We do not guarantee to throw either `BoundsError` or `BoundsErrorLight` when accessing an out-of-bound index of a `FixedSizeArray`.
    The exact exception raised may change at any point without breaking the semantic versioning compatibility contract.

!!! note

    `BoundsErrorLight` is an alias for `LightBoundsError` from the [LightBoundsErrors.jl](https://github.com/JuliaArrays/LightBoundsErrors.jl) package.
