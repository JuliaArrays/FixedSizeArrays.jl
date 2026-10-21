# Release notes

The format is based on [Keep a Changelog](https://keepachangelog.com/en/1.0.0/), and this project adheres to [Semantic Versioning](https://semver.org/spec/v2.0.0.html).

## Version [v1.5.0](https://github.com/JuliaArrays/FixedSizeArrays.jl/releases/tag/v1.5.0) - Unreleased

### Added

* Static dimensions: `FixedSizeArray{T,N,Mem}` is now an alias for `FixedSizeArrayDims{T,N,NTuple{N,Int},Mem}`, where the new public (but not exported) type `FixedSizeArrays.FixedSizeArrayDims{T,N,D,Mem}` is an `N`-dimensional array over a memory backend `Mem<:DenseVector{T}` whose dimensions are a tuple of type `D` of `Int`s, `Val`s and `StaticInt`s from [Static.jl](https://github.com/SciML/Static.jl) (e.g. `(Val(5), Val(4), n)` or `(static(5), static(4), n)`), so that any dimension can be a compile-time constant while the others are set at runtime. Linear indices, strides and sizes are computed at compile time as far as the dimensions are static. An array with static dimensions can be constructed with `FixedSizeArrays.wrap(FixedSizeArrayDims, parent, dims)`, with `collect_as`, or with the constructors of `FixedSizeArrayDims` (`undef` or converting from another array); the dimensions given to a `FixedSizeArray` constructor are converted to `Int`s. It supports the strided arrays interface, `similar`, `copy`, `fill!`, `copyto!`, broadcasting (returning `FixedSizeArray`s), and the package extensions with MPI.jl, Adapt.jl (adapting the memory backend for any memory backend), and the new package extension with [StaticArrayInterface.jl](https://github.com/JuliaArrays/StaticArrayInterface.jl), which provides static sizes and strides to e.g. StrideArraysCore.jl, LoopVectorization.jl and Polyester.jl. `Val` dimensions need no other package, while `StaticInt` dimensions are available with the new package extension with Static.jl, which is a weak dependency.
* `FixedSizeArrays.wrap` now accepts at most one dimension given as `:`, which is computed from the length of the parent.
* With StaticArrayInterface.jl, `StaticArrayInterface.device` of a `FixedSizeArray` is now the one of its memory backend, e.g. a GPU memory backend is not treated as a CPU pointer by LoopVectorization.jl.

### Changed

* `FixedSizeArray{T,N,Mem}` is an alias for `FixedSizeArrayDims{T,N,NTuple{N,Int},Mem}` (see above). It is still a concrete type, with the same layout and methods as before, but `FixedSizeArray` itself is no longer the name of a `struct`: `Base.typename(FixedSizeArray)` is `FixedSizeArrayDims`, and `FixedSizeArray{T,N}` with an invalid `N` (e.g. `FixedSizeArray{Int,-1}`) now throws an error when the type is formed.
* Cartesian `getindex` and `setindex!` with one `Int` index per dimension compute the linear index in Horner form, `i1 + n1 * ((i2 - 1) + n2 * (...))`, instead of a sum of strided terms, `1 + (i1 - 1) + n1 * (i2 - 1) + (n1 * n2) * (...)`. This needs fewer integer operations, e.g. in kernels accessing all elements `u[v, i, j, k]` along the first dimension. ([#217](https://github.com/JuliaArrays/FixedSizeArrays.jl/pull/217))

### Fixed

* Improved performance of loops over `eachindex` which store into a `FixedSizeArray` on Julia v1.13 and later: `eachindex` of a `FixedSizeArray` no longer clamps its length to non-negative values, which lets LLVM elide the bounds checks of these loops and vectorize them. ([#200](https://github.com/JuliaArrays/FixedSizeArrays.jl/issues/200))

## Version [v1.4.2](https://github.com/JuliaArrays/FixedSizeArrays.jl/releases/tag/v1.4.2) - 2026-10-13

### Added

* No public-facing changes in this version, only new FileCheck-based tests added ([#204](https://github.com/JuliaArrays/FixedSizeArrays.jl/pull/204)).

## Version [v1.4.1](https://github.com/JuliaArrays/FixedSizeArrays.jl/releases/tag/v1.4.1) - 2026-10-13

### Fixed

* Improved performance on Julia v1.13 and later of loops after a check of the size of a `FixedSizeArray` argument, such as `axes(a) == axes(b) || throw(...)`: `size` now loads all dimensions at once, so that LLVM can use the checked dimensions in the loop instead of reloading them. ([#201](https://github.com/JuliaArrays/FixedSizeArrays.jl/issues/201), [#202](https://github.com/JuliaArrays/FixedSizeArrays.jl/pull/202))

## Version [v1.4.0](https://github.com/JuliaArrays/FixedSizeArrays.jl/releases/tag/v1.4.0) - 2026-10-13

### Added

* New public function `FixedSizeArrays.wrap` to construct a `FixedSizeArray` around an existing parent `DenseVector` (e.g. a `Memory` or a `Vector`) without copying it. ([#196](https://github.com/JuliaArrays/FixedSizeArrays.jl/pull/196))
* Package extension with [`MPI.jl`](https://github.com/JuliaParallel/MPI.jl): `FixedSizeArray`s (and `SubArray`s of them) can now be used directly as communication buffers in MPI operations, e.g. `MPI.Send`, `MPI.Recv!`, `MPI.Bcast!`, `MPI.Allreduce!` and `MPI.Gatherv!`, by defining `MPI.Buffer(::FixedSizeArray)` and the conversion to `MPI.MPIPtr`. Requires MPI.jl v0.20.28 or later. ([#198](https://github.com/JuliaArrays/FixedSizeArrays.jl/pull/198))

### Changed

* `BoundsErrorLight` is now an alias for `LightBoundsError` from the [`LightBoundsErrors.jl`](https://github.com/JuliaArrays/LightBoundsErrors.jl) package, instead of being a type defined in this package. As a consequence, the properties of the exception object and the error message have changed. ([#176](https://github.com/JuliaArrays/FixedSizeArrays.jl/pull/176))
* All out-of-bounds accesses to a `FixedSizeArray` going through `checkbounds` now throw `LightBoundsError`, not only those indexing with a single integer; the other cases used to throw `BoundsError`. ([#175](https://github.com/JuliaArrays/FixedSizeArrays.jl/issues/175), [#176](https://github.com/JuliaArrays/FixedSizeArrays.jl/pull/176))
* The validation of the dimensions when constructing a `FixedSizeArray` is now implemented with the [`CheckedSizeProduct.jl`](https://github.com/JuliaArrays/CheckedSizeProduct.jl) package. Invalid dimensions still throw an `ArgumentError`, but with a single error message which includes the offending shape. ([#89](https://github.com/JuliaArrays/FixedSizeArrays.jl/issues/89), [#90](https://github.com/JuliaArrays/FixedSizeArrays.jl/pull/90))

### Fixed

* `isassigned(::FixedSizeArray, ::Int)` now checks the bounds on the `FixedSizeArray` itself, instead of relying on the parent to do it, giving correct results also with underlying storage types which do not check bounds. ([#188](https://github.com/JuliaArrays/FixedSizeArrays.jl/pull/188))
* `copy!(::FixedSizeVector, ::AbstractVector)` now throws a descriptive `ArgumentError` when the axes of the two arrays do not match, instead of a `MethodError` caused by trying to call `resize!` on the destination. ([#194](https://github.com/JuliaArrays/FixedSizeArrays.jl/issues/194), [#195](https://github.com/JuliaArrays/FixedSizeArrays.jl/pull/195))

## Version [v1.3.0](https://github.com/JuliaArrays/FixedSizeArrays.jl/releases/tag/v1.3.0) - 2026-02-14

### Changed

* The package extensions were renamed from `RandomExt` and `AdaptExt` to `FixedSizeArraysRandomExt` and `FixedSizeArraysAdaptExt`, respectively, to avoid clashes with similarly named extensions in other packages. ([#162](https://github.com/JuliaArrays/FixedSizeArrays.jl/pull/162))
* Improved performance of iterating over a `FixedSizeArray`: constant propagation through `iterate` now allows the compiler to unroll loops over arrays with statically known size. ([#158](https://github.com/JuliaArrays/FixedSizeArrays.jl/pull/158))
* Improved performance of destructuring a `FixedSizeArray` with a slurping pattern (e.g. `a, rest... = array`), by adding a method to `Base.rest`. ([#166](https://github.com/JuliaArrays/FixedSizeArrays.jl/pull/166))

### Fixed

* Indexing a multidimensional `FixedSizeArray` with a single `Colon` (`A[:]`) now always returns a vector, instead of an array with the same shape as `A`. ([#161](https://github.com/JuliaArrays/FixedSizeArrays.jl/pull/161))
* Fixed type inference of `collect_as` in Julia v1.12 and later. ([#160](https://github.com/JuliaArrays/FixedSizeArrays.jl/issues/160), [#172](https://github.com/JuliaArrays/FixedSizeArrays.jl/pull/172))

## Version [v1.2.0](https://github.com/JuliaArrays/FixedSizeArrays.jl/releases/tag/v1.2.0) - 2025-08-07

### Added

* Pkgextension with [`Adapt.jl`](https://github.com/JuliaGPU/Adapt.jl), for converting other arrays to `FixedSizeArray`. ([#149](https://github.com/JuliaArrays/FixedSizeArrays.jl/pull/149))

## Version [v1.1.0](https://github.com/JuliaArrays/FixedSizeArrays.jl/releases/tag/v1.1.0) - 2025-08-05

### Added

* Support for broadcasting of view of `FixedSizeArray`s. ([#146](https://github.com/JuliaArrays/FixedSizeArrays.jl/issues/146), [#147](https://github.com/JuliaArrays/FixedSizeArrays.jl/pull/147))

## Version [v1.0.1](https://github.com/JuliaArrays/FixedSizeArrays.jl/releases/tag/v1.0.1) - 2025-06-08

### Fixed

* `rand!(A::FixedSizeArray)` now returns `A`, instead of its parent. ([#139](https://github.com/JuliaArrays/FixedSizeArrays.jl/issues/139), [#140](https://github.com/JuliaArrays/FixedSizeArrays.jl/pull/140))

## Version [v1.0.0](https://github.com/JuliaArrays/FixedSizeArrays.jl/releases/tag/v1.0.0) - 2025-06-04

Initial release.
