## `FixedSizeArrays.jl`

[![docs stable](https://img.shields.io/badge/docs-stable-blue.svg)](https://juliaarrays.github.io/FixedSizeArrays.jl/stable)
[![docs dev](https://img.shields.io/badge/docs-dev-blue.svg)](https://juliaarrays.github.io/FixedSizeArrays.jl/dev)
[![Aqua QA](https://raw.githubusercontent.com/JuliaTesting/Aqua.jl/master/badge.svg)](https://github.com/JuliaTesting/Aqua.jl)
[![Build Status](https://github.com/JuliaArrays/FixedSizeArrays.jl/actions/workflows/UnitTests.yml/badge.svg?branch=main)](https://github.com/JuliaArrays/FixedSizeArrays.jl/actions/workflows/UnitTests.yml?query=branch%3Amain)
[![Coverage](https://codecov.io/gh/JuliaArrays/FixedSizeArrays.jl/branch/main/graph/badge.svg)](https://codecov.io/gh/JuliaArrays/FixedSizeArrays.jl)
[![PkgEval](https://JuliaCI.github.io/NanosoldierReports/pkgeval_badges/F/FixedSizeArrays.svg)](https://JuliaCI.github.io/NanosoldierReports/pkgeval_badges/F/FixedSizeArrays.html)

[`FixedSizeArrays.jl`](https://github.com/JuliaArrays/FixedSizeArrays.jl) is a package for the [Julia programming language](https://julialang.org/) which implements mutable fixed-size arrays: arrays whose size is set once, at construction, and cannot change afterwards.
The size can be fixed at two different stages:

* at runtime, with `FixedSizeArray`: the size is a constant field of the array, so arrays of all sizes share the same type and compiled code, and the compiler can [constant-propagate](https://en.wikipedia.org/wiki/Constant_folding) it when it is known at compile time;
* at compile time, dimension by dimension, with `FixedSizeArrays.FixedSizeArrayDims`: each dimension can be static, a compile-time constant which is part of the type, or dynamic, an `Int` stored in the array.
  Static dimensions require the [Static.jl](https://github.com/SciML/Static.jl) package, a weak dependency.

`FixedSizeArray` is an alias for `FixedSizeArrayDims` with all dimensions dynamic, so everything which works for a `FixedSizeArray` also works for an array with static dimensions.

FixedSizeArrays supports the [standard array interfaces](https://docs.julialang.org/en/v1/manual/interfaces/#man-interface-array), so things like broadcasting, matrix multiplication, other linear algebra operations, `similar`, `copyto!` or `map` should just work.

Find more in the manual:

* [documentation of latest stable version](https://juliaarrays.github.io/FixedSizeArrays.jl/stable)
* [documentation of in-development version](https://juliaarrays.github.io/FixedSizeArrays.jl/dev)

## Videos about `FixedSizeArrays.jl`

* [FixedSizeArrays.jl: What Array probably should have been](https://www.youtube.com/watch?v=xWo-ttbSHgw), talk at JuliaCon '25 by Mosè Giordano, Oscar Smith, and Neven Sajko.
