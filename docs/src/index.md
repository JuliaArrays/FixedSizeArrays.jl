# FixedSizeArrays.jl

[`FixedSizeArrays.jl`](https://github.com/JuliaArrays/FixedSizeArrays.jl) is a package for the [Julia programming language](https://julialang.org/) which implements mutable fixed-size arrays ([`FixedSizeArray`](@ref)), which means the length of the array can not change after construction and is amenable to be [constant-propagated](https://en.wikipedia.org/wiki/Constant_folding) at compile-time when possible.

FixedSizeArrays supports the [standard array interfaces](https://docs.julialang.org/en/v1/manual/interfaces/#man-interface-array), so things like broadcasting, matrix multiplication, other linear algebra operations, `similar`, `copyto!` or `map` should just work.

Use the constructors to convert from other array types.
Use `collect_as` from the [Collects.jl](https://github.com/JuliaCollections/Collects.jl) package to convert from arbitrary iterators.

## Comparison with other array types

### `Array` from `Base`

While `Base`'s `Array` is a convenient multi-dimensional container for continuous data, a fixed-size array is what is employed in most linear algebra applications, where the size of [tensors](https://en.wikipedia.org/wiki/Tensor) typically does not change.
Here is a comparison between `Array` and `FixedSizeArray`:

* `Array` is [dynamic](https://en.wikipedia.org/wiki/Dynamic_array) (or growable) which means the size of an `Array` can be changed after construction.  The size of a `FixedSizeArray` cannot be changed;
* the above property of `FixedSizeArray`s enables certain compiler optimizations in some cases where the size of a `FixedSizeArray` is known at compile-time.
  The dynamic property of `Array`s limits the optimization opportunities, which makes `Array`s less likely to be optimized as efficiently as `FixedSizeArray`s.

You can see a showcase of performance/memory advantages of `FixedSizeArray` over `Array` in [this evolving discussion](https://github.com/JuliaArrays/FixedSizeArrays.jl/discussions/62).

### `MArray` from `StaticArrays.jl`

`FixedSizeArray` is an alternative implementation to [`MArray`](https://juliaarrays.github.io/StaticArrays.jl/stable/api/#StaticArraysCore.MArray) from [`StaticArrays.jl`](https://github.com/JuliaArrays/StaticArrays.jl).

Main differences between `FixedSizeArray` and `MArray` are:

* `FixedSizeArray` is based on the `Memory` type introduced in Julia v1.11, `MArray` is backed by tuples, which are not efficient for long arrays (``\gg 100`` elements);
* the size of the array is part of the type parameters of `MArray`, this isn't the case for `FixedSizeArray`, where the size is only a constant field of the data structure.

There are applications where it may be desirable to dispatch on the size of an array, in that case `FixedSizeArray` would not be a good fit and you may consider using `MArray` instead.

`FixedSizeArray` is not a direct replacement for `SArray`: in addition to having the size as part of the type parameters, `SArray` is immutable, instead the elements of a `FixedSizeArray` can be changed at any point.

## What about the other package with same name?

This package is not related to a package with the same name by [`@SimonDanisch`](https://github.com/SimonDanisch/FixedSizeArrays.jl).
That, earlier, package was one of the StaticArrays.jl-like packages in the pre-v1 days of Julia: <https://github.com/SimonDanisch/FixedSizeArrays.jl>.
