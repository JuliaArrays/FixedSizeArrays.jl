# Usage

This page covers the two kinds of arrays of the package in turn: [`FixedSizeArray`](@ref), whose size is set at runtime, and then [`FixedSizeArrays.FixedSizeArrayDims`](@ref), the underlying type in which some dimensions can be set at compile time.
Everything in the first section applies to both, since a `FixedSizeArray` is a `FixedSizeArrayDims` with all dimensions dynamic.

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
2×2 FixedSizeArrays.FixedSizeArrayDims{Int64, 2, Tuple{Int64, Int64}, Memory{Int64}}:
 10  20
 30  14

julia> FixedSizeArray{Float64}(arr)  # construct from an `AbstractArray` value while converting element type
2×2 FixedSizeArrays.FixedSizeArrayDims{Float64, 2, Tuple{Int64, Int64}, Memory{Float64}}:
 10.0  20.0
 30.0  14.0
```

!!! note "The printed type"

    Depending on the version of Julia, the type of a `FixedSizeArray` may be printed as `FixedSizeArrays.FixedSizeArrayDims{Int64, 2, Tuple{Int64, Int64}, Memory{Int64}}`, as above, rather than as `FixedSizeMatrix{Int64, Memory{Int64}}`.
    `FixedSizeArray{T,N,Mem}` is an alias for `FixedSizeArrayDims{T,N,NTuple{N,Int},Mem}`: the third parameter is the type of the dimensions of the array, here two `Int`s, i.e. two dimensions set at runtime.
    See [Static dimensions: the `FixedSizeArrayDims` type](@ref) for the arrays in which some of them are `StaticInt`s instead.

To "convert" a wrapper type like `Adjoint` while preserving the wrapper structure, you can use `Adapt`:

```jldoctest
julia> using FixedSizeArrays, Adapt

julia> arr = transpose([10 20; 30 14])
2×2 transpose(::Matrix{Int64}) with eltype Int64:
 10  30
 20  14

julia> adapt(FixedSizeArray, arr)
2×2 transpose(::FixedSizeArrays.FixedSizeArrayDims{Int64, 2, Tuple{Int64, Int64}, Memory{Int64}}) with eltype Int64:
 10  30
 20  14
```

For 1- and 2-dimensional arrays you can use the aliases [`FixedSizeVector`](@ref) and [`FixedSizeMatrix`](@ref), respectively:

```@repl
using FixedSizeArrays
v = FixedSizeVector{Float64}(undef, 2)
M = FixedSizeMatrix{Float64}(undef, 3, 3)
```

The [array literals syntax](https://docs.julialang.org/en/v1/manual/arrays/#man-array-literals) `[A, B, C, ...]` is limited to `Base`'s `Array`: to construct a `FixedSizeArray` from any iterator, use `collect_as` from Collects.jl, see [The `collect_as` utility](@ref).

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

The memory backend of an array is returned by `parent`.

### Wrapping an existing memory backend

The constructors above copy the elements into a new memory backend.
[`FixedSizeArrays.wrap`](@ref) instead returns a `FixedSizeArray` which uses an existing `DenseVector`, e.g. a `Memory` or a `Vector`, as its memory backend, without copying:

```jldoctest
julia> using FixedSizeArrays

julia> v = collect(1:6);

julia> a = FixedSizeArrays.wrap(FixedSizeArray, v, (2, 3));

julia> size(a), parent(a) === v
((2, 3), true)

julia> a[2, 3] = 60;  # writes to `v`

julia> v[6]
60
```

The length of the parent must be equal to the product of the dimensions.
At most one dimension can be given as `:`, which is then computed from the length of the parent.
Changing the length of the parent while the wrapped array is in use, e.g. with `resize!` when the parent is a `Vector`, is not allowed: wrap the parent again after resizing it instead.

## Static dimensions: the `FixedSizeArrayDims` type

### Dynamic and static dimensions

Many arrays in numerical codes have a shape that is partly known in advance.
Take the solution of a discontinuous Galerkin method in three dimensions, stored as `u[variable, i, j, k, element]`: the number of variables (five for the compressible Euler equations) and of nodes per direction (say four) are fixed by the choice of the equations and the polynomial degree, but the number of elements depends on the mesh, and may change while the simulation runs, e.g. with adaptive mesh refinement.
With a plain `Array` or `FixedSizeArray` the compiler knows none of these sizes, so every index computation multiplies with sizes loaded from memory, and loops over the small dimensions cannot be fully unrolled.

For such arrays, this package lets you choose, for each dimension, whether it is *dynamic* or *static*:

* a dynamic dimension is an `Int` stored in the array, as all dimensions of a `FixedSizeArray`;
* a static dimension is a `StaticInt` from [Static.jl](https://github.com/SciML/Static.jl): it is part of the type of the array and takes no space in it, so the compiler can compute the size in that dimension, and the strides and linear indices involving it, at compile time.

The dimensions are the third parameter of the underlying type [`FixedSizeArrays.FixedSizeArrayDims{T,N,D,Mem}`](@ref FixedSizeArrays.FixedSizeArrayDims): `D` is a tuple type of `N` `Int`s and `StaticInt`s, in any order.
`FixedSizeArray{T,N,Mem}` is the alias `FixedSizeArrayDims{T,N,NTuple{N,Int},Mem}`, so all `FixedSizeArray`s are `FixedSizeArrayDims`s, but the arrays with static dimensions are not `FixedSizeArray`s.
`FixedSizeArrayDims` is public, but not exported: use it as `FixedSizeArrays.FixedSizeArrayDims`, or import it with `using FixedSizeArrays: FixedSizeArrayDims`.
Static.jl is a weak dependency of FixedSizeArrays.jl, so static dimensions are only available once you load it yourself, e.g. with `using Static: static`.

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

Each combination of static dimensions is a different type, so the code using the array is compiled for each one: make static only the dimensions which are known in advance and matter for performance, such as the small ones in the example above.

### Constructing arrays with static dimensions

The static dimensions cannot be derived from an array or an iterator, so an array with static dimensions is constructed either with dimensions given as `StaticInt`s, or with a type which specifies the dimensions type `D`:

* the `undef` constructor `FixedSizeArrayDims{T}(undef, dims)` makes static the dimensions given as `static(n)`, and dynamic those given as `Int`s, as in the example above;
* the converting constructor and [`collect_as`](@ref "The `collect_as` utility") need the dimensions type: `FixedSizeArrayDims{Float64, 2, Tuple{StaticInt{2}, Int}}(arr)` copies the matrix `arr` into an array whose first dimension is static, and `collect_as(FixedSizeArrayDims{Int, 1, Tuple{StaticInt{3}}}, iter)` collects three elements into a vector of static length;
* [`FixedSizeArrays.wrap`](@ref) wraps an existing memory backend without copying, with the dimensions given as `StaticInt`s.

```jldoctest
julia> using FixedSizeArrays: FixedSizeArrays, FixedSizeArrayDims

julia> using Static: static, StaticInt

julia> using Collects: collect_as

julia> a = FixedSizeArrayDims{Float64, 2, Tuple{StaticInt{2}, Int}}([1 2 3; 4 5 6]);  # copy from a matrix

julia> size(a), size(a, 1)
((2, 3), 2)

julia> b = collect_as(FixedSizeArrayDims{Int, 1, Tuple{StaticInt{3}}}, (i^2 for i in 1:3));  # from an iterator

julia> b == [1, 4, 9]
true

julia> v = collect(1.0:12.0);

julia> c = FixedSizeArrays.wrap(FixedSizeArrayDims, v, (static(2), static(3), :));  # wrap, computing the last dimension

julia> size(c), parent(c) === v
((2, 3, 2), true)

julia> d = FixedSizeArrayDims{Int}(undef, (5, static(3)));  # a static dimension after a dynamic one

julia> size(d)
(5, 3)
```

When the type of the array specifies the dimensions type, the dimensions given to the constructor or to `wrap` are converted to it: `FixedSizeArray{Float64}(undef, (static(2), 3))` has two dynamic dimensions, since a `FixedSizeArray` has only dynamic dimensions, while `FixedSizeArrayDims{Float64}(undef, (static(2), 3))` has a static one.

### Operations which preserve the static dimensions

Operations which return an array with the same shape as the input keep its static dimensions, so the result has the same type: `similar(a)`, `similar(a, T)`, `similar(typeof(a), size)` and `copy`, as well as the in-place operations such as `copyto!`, `fill!`, in-place broadcasting with `.=` or `rand!`.
Operations whose result has a shape computed at runtime return a `FixedSizeArray` with all dimensions dynamic, so that the type of their result does not depend on the value of the shape: `similar(a, T, size)`, out-of-place broadcasting, `map`, reductions, `reshape`, `vec` and indexing with ranges or colons.

```jldoctest
julia> using FixedSizeArrays

julia> using FixedSizeArrays: FixedSizeArrayDims

julia> using Static: static

julia> a = FixedSizeArrayDims{Float64}(undef, (static(2), 3));

julia> typeof(similar(a)) === typeof(a)  # keeps the static dimension
true

julia> typeof(a .+ 1) === FixedSizeArray{Float64, 2, Memory{Float64}}  # all dimensions dynamic
true

julia> typeof(a[:, 1]) === FixedSizeVector{Float64, Memory{Float64}}
true
```

The static dimensions are also available to other packages through the package extensions for Adapt.jl and StaticArrayInterface.jl, see [Package extensions](@ref).

## Package extensions

The package has [extensions](https://pkgdocs.julialang.org/v1/creating-packages/#Conditional-loading-of-code-in-packages-(Extensions)) which are loaded automatically together with the following packages:

* [Adapt.jl](https://github.com/JuliaGPU/Adapt.jl): `adapt(FixedSizeArray, x)` converts the arrays in `x` to `FixedSizeArray`s, also inside wrappers such as `Adjoint`, see the [constructors](@ref "Constructors") above.
  Conversely, `adapt(to, a)` converts a `FixedSizeArray` `a` as a whole, e.g. `adapt(CuArray, a)` is a `CuArray`, while for an array with static dimensions it adapts the memory backend and wraps it again, so that `adapt(CuArray, a)` is a `FixedSizeArrayDims` backed by a `CuVector`, with the same static dimensions.
  `FixedSizeArrayDims` itself is not an `adapt` target: `adapt(FixedSizeArrayDims, x)` returns `x` unchanged, since the static dimensions can not be derived from `x`, use the [constructors](@ref "Constructing arrays with static dimensions") instead;
* [MPI.jl](https://github.com/JuliaParallel/MPI.jl) (v0.20.28 or later): the arrays of this package, and `SubArray`s of them, can be used directly as communication buffers in MPI operations, e.g. `MPI.Send`, `MPI.Recv!`, `MPI.Bcast!`, `MPI.Allreduce!` and `MPI.Gatherv!`;
* [Random.jl](https://docs.julialang.org/en/v1/stdlib/Random/): `rand!(a)` fills the array in place;
* [Static.jl](https://github.com/SciML/Static.jl): enables the [static dimensions](@ref "Static dimensions: the `FixedSizeArrayDims` type");
* [StaticArrayInterface.jl](https://github.com/JuliaArrays/StaticArrayInterface.jl), together with Static.jl: the static sizes and strides of an array with static dimensions, and the device of its memory backend, are made available to the packages which use this interface, such as LoopVectorization.jl, StrideArraysCore.jl and Polyester.jl.

## The `collect_as` utility

The [array literals syntax](https://docs.julialang.org/en/v1/manual/arrays/#man-array-literals) `[A, B, C, ...]` to construct arrays is limited to `Base`'s `Array` and cannot be extended to custom array types.
The package [Collects.jl](https://github.com/JuliaCollections/Collects.jl) provides a convenient function, `collect_as`, to overcome this limitation and construct `FixedSizeArray`s out of any iterable:

```jldoctest
julia> iter = (i for i ∈ 7:9 if i≠8);

julia> using FixedSizeArrays, Collects

julia> collect_as(FixedSizeArray, iter)  # construct from an arbitrary iterator
2-element FixedSizeArrays.FixedSizeArrayDims{Int64, 1, Tuple{Int64}, Memory{Int64}}:
 7
 9

julia> collect_as(FixedSizeArray{Float64}, iter)  # construct from an arbitrary iterator while converting element type
2-element FixedSizeArrays.FixedSizeArrayDims{Float64, 1, Tuple{Int64}, Memory{Float64}}:
 7.0
 9.0

julia> collect_as(FixedSizeVectorDefault, (3.14, -4.2, 2.68))  # construct from a tuple
3-element FixedSizeArrays.FixedSizeArrayDims{Float64, 1, Tuple{Int64}, Memory{Float64}}:
  3.14
 -4.2
  2.68
```

`collect_as` also constructs arrays with static dimensions, given a `FixedSizeArrayDims` type with the dimensions type, see [Constructing arrays with static dimensions](@ref).
See the Collects.jl Readme for more information.

## `BoundsErrorLight` exception

To facilitate the [escape analysis](https://en.wikipedia.org/wiki/Escape_analysis) of the arrays of this package, accessing an out-of-bound index of these arrays raises a `BoundsErrorLight` exception when possible.
This exception type does not store the entire array for reporting the error message, thus enabling more performance optimizations compared to arrays which throw [`BoundsError`](https://docs.julialang.org/en/v1/base/base/#Core.BoundsError) exceptions.

!!! warning "Compatibility of exception raised"

    We do not guarantee to throw either `BoundsError` or `BoundsErrorLight` when accessing an out-of-bound index of a `FixedSizeArray`.
    The exact exception raised may change at any point without breaking the semantic versioning compatibility contract.

!!! note

    `BoundsErrorLight` is an alias for `LightBoundsError` from the [LightBoundsErrors.jl](https://github.com/JuliaArrays/LightBoundsErrors.jl) package.
