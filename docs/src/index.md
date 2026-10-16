# FixedSizeArrays.jl

[`FixedSizeArrays.jl`](https://github.com/JuliaArrays/FixedSizeArrays.jl) is a package for the [Julia programming language](https://julialang.org/) which implements mutable fixed-size arrays: arrays whose size is set once, at construction, and cannot change afterwards.
The size can be fixed at two different stages:

* at runtime, with [`FixedSizeArray`](@ref): the size is a constant field of the array, so arrays of all sizes have the same type and share the same compiled code.
  Since the size cannot change, the compiler is free to [constant-propagate](https://en.wikipedia.org/wiki/Constant_folding) it whenever it is known at compile time;
* at compile time, dimension by dimension, with [`FixedSizeArrays.FixedSizeArrayDims`](@ref): each dimension can be *static*, a compile-time constant which is part of the type and takes no space in the array, or *dynamic*, an `Int` stored in the array as in a `FixedSizeArray`.
  The compiler then computes the sizes, strides and linear indices along the static dimensions at compile time, in any function which receives the array.

These are not two separate array types: `FixedSizeArray` is an alias for `FixedSizeArrayDims` with all dimensions dynamic, so everything which works for a `FixedSizeArray` also works for an array with static dimensions.
Static dimensions are opt-in: they require the [Static.jl](https://github.com/SciML/Static.jl) package, which is a weak dependency of FixedSizeArrays.jl.

Both kinds of arrays support the [standard array interfaces](https://docs.julialang.org/en/v1/manual/interfaces/#man-interface-array), so things like broadcasting, matrix multiplication, other linear algebra operations, `similar`, `copyto!` or `map` should just work.
The elements are stored contiguously in a *memory backend*, which can be any `DenseVector`, e.g. a GPU vector.
The default memory backend is the [`Memory`](https://docs.julialang.org/en/v1/base/arrays/#Core.Memory) type introduced in Julia v1.11, the same type which backs `Base`'s `Array`.

```@repl
using FixedSizeArrays
a = FixedSizeArray([1.0 2 3; 4 5 6]) # size fixed at runtime: (2, 3)
a[1, 1] = 10
using Static: static
using FixedSizeArrays: FixedSizeArrayDims
b = FixedSizeArrayDims{Float64}(undef, (static(2), 3)); # the first dimension is part of the type
copyto!(b, a)
size(b, 1) # a compile-time constant
```

## Which one should I use?

Start with `FixedSizeArray`.
It is a replacement for `Array` wherever you do not need to resize an array, which is the case for most numerical code, with the same contiguous storage but more optimization opportunities for the compiler, see the [comparison with `Array`](@ref "`Array` from `Base`") below.
The [Usage](@ref "The `FixedSizeArray` type") page shows how to construct it from scratch, from other arrays, with [`FixedSizeArrays.wrap`](@ref) around an existing memory backend, or with `collect_as` from any iterator.

Reach for the static dimensions of `FixedSizeArrayDims` when some dimensions of an array are small and known in advance, e.g. the number of components of a vector field, and the loops over them should be unrolled or the index computations along them hoisted at compile time, or when you need to dispatch on the size of a dimension.
The price is that each combination of static dimensions is a different type, so the code using the array is compiled for each one, as with StaticArrays.jl.
The remaining dimensions, e.g. the large one which depends on the input data, stay dynamic.
See [Static dimensions: the `FixedSizeArrayDims` type](@ref) for how to construct these arrays and which operations preserve the static dimensions.

## Comparison with other array types

Here is a quick summary of high-level differences between some different array types:

|                           | Size set at...                         | Data backend        | Growable | Mutable elements      |
|---------------------------|----------------------------------------|---------------------|----------|-----------------------|
| `Base.Array`              | runtime                                | `Memory`            | ✅        | ✅                     |
| `FixedSizeArray`          | runtime                                | any `DenseVector`   | ❌        | ✅                     |
| `FixedSizeArrayDims`      | compile time or runtime, per dimension | any `DenseVector`   | ❌        | ✅                     |
| `MArray`                  | compile time                           | `Tuple`             | ❌        | ✅                     |
| `SArray`                  | compile time                           | `Tuple`             | ❌        | ❌                     |
| `HybridArray`             | compile time or runtime, per dimension | any `AbstractArray` | ❌        | depends on the parent |
| `StrideArray`, `PtrArray` | compile time or runtime, per dimension | a pointer           | ❌        | ✅                     |

In the following sections we elaborate on the properties and strengths of the various types.

### `Array` from `Base`

While `Base`'s `Array` is a convenient multi-dimensional container for contiguous data, a fixed-size array is what is employed in most linear algebra applications, where the size of [tensors](https://en.wikipedia.org/wiki/Tensor) typically does not change.
Here is a comparison between `Array` and the arrays of this package:

* `Array` is [dynamic](https://en.wikipedia.org/wiki/Dynamic_array) (or growable) which means the size of an `Array` can be changed after construction.
  The size of a `FixedSizeArray` or a `FixedSizeArrayDims` cannot be changed;
* the above property enables certain compiler optimizations in some cases where the size of a `FixedSizeArray` is known at compile time, e.g. by constant propagation.
  The dynamic property of `Array`s limits the optimization opportunities, which makes `Array`s less likely to be optimized as efficiently as `FixedSizeArray`s;
* with a `FixedSizeArrayDims`, the static dimensions do not depend on constant propagation: they are part of the type, so the compiler knows them in any function which receives the array.

You can see a showcase of performance/memory advantages of `FixedSizeArray` over `Array` in [this evolving discussion](https://github.com/JuliaArrays/FixedSizeArrays.jl/discussions/62).

### `MArray` and `SArray` from `StaticArrays.jl`

`FixedSizeArray` is an alternative implementation to [`MArray`](https://juliaarrays.github.io/StaticArrays.jl/stable/api/#StaticArraysCore.MArray) from [`StaticArrays.jl`](https://github.com/JuliaArrays/StaticArrays.jl).

Main differences between the arrays of this package and `MArray` are:

* `FixedSizeArray` and `FixedSizeArrayDims` are backed by a `DenseVector`, by default the `Memory` type introduced in Julia v1.11, while `MArray` is backed by tuples, which are not efficient for long arrays (``\gg 100`` elements);
* the whole size of an `MArray` is part of its type parameters, while the size of a `FixedSizeArray` is only a constant field of the data structure, so arrays of all sizes have the same type and share the same compiled code.
  `FixedSizeArrayDims` sits in between: you choose which dimensions are part of the type, e.g. the small dimensions known in advance, while the large ones are set at runtime.

There are applications where it may be desirable to dispatch on the size of an array.
That is not possible with `FixedSizeArray`, while it is with the static dimensions of a `FixedSizeArrayDims`, or with `MArray`.
However, StaticArrays.jl implements many operations, e.g. linear algebra with small matrices, with fully unrolled code whose results have static sizes too, while the results of out-of-place operations on a `FixedSizeArrayDims`, such as broadcasting, `map` or reductions, are `FixedSizeArray`s with all dimensions dynamic.
For small arrays whose size should propagate through all operations, `MArray` and `SArray` are a better fit.

Neither `FixedSizeArray` nor `FixedSizeArrayDims` is a direct replacement for `SArray`: in addition to having the size as part of the type parameters, `SArray` is immutable, instead the elements of the arrays of this package can be changed at any point.

### Other arrays with static dimensions: `HybridArray` and `StrideArray`

Other packages also mix static and dynamic dimensions, as `FixedSizeArrayDims` does:

* [HybridArrays.jl](https://github.com/JuliaArrays/HybridArrays.jl) builds on StaticArrays.jl, using `StaticArrays.Dynamic()` for the dynamic dimensions, so that e.g. indexing along the static dimensions can return static arrays.
  It wraps any parent array, so it can be non-contiguous, while a `FixedSizeArrayDims` is always a `DenseArray` backed by the whole memory backend.
* [StrideArraysCore.jl](https://github.com/JuliaSIMD/StrideArraysCore.jl) uses the same `Int`/`StaticInt` tuples for the sizes, and adds static strides and offsets, for LoopVectorization.jl.
  Its arrays are built around a pointer (a `PtrArray` must be kept alive with `GC.@preserve`), while a `FixedSizeArrayDims` holds a reference to its memory backend, so the garbage collector knows about it, and the memory backend can be a GPU vector.
  The [package extension for StaticArrayInterface.jl](@ref "Package extensions") makes `FixedSizeArrayDims`s work with StrideArraysCore.jl and LoopVectorization.jl, too.

## What about the other package with same name?

This package is not related to a package with the same name by [`@SimonDanisch`](https://github.com/SimonDanisch/FixedSizeArrays.jl).
That, earlier, package was one of the StaticArrays.jl-like packages in the pre-v1 days of Julia: <https://github.com/SimonDanisch/FixedSizeArrays.jl>.

## Videos about `FixedSizeArrays.jl`

* [FixedSizeArrays.jl: What Array probably should have been](https://www.youtube.com/watch?v=xWo-ttbSHgw), talk at JuliaCon '25 by Mosè Giordano, Oscar Smith, and Neven Sajko.
