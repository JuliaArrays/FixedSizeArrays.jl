## `FixedSizeArrays.jl`

[![docs stable](https://img.shields.io/badge/docs-stable-blue.svg)](https://juliaarrays.github.io/FixedSizeArrays.jl/stable)
[![docs dev](https://img.shields.io/badge/docs-dev-blue.svg)](https://juliaarrays.github.io/FixedSizeArrays.jl/dev)
[![Aqua QA](https://raw.githubusercontent.com/JuliaTesting/Aqua.jl/master/badge.svg)](https://github.com/JuliaTesting/Aqua.jl)
[![Build Status](https://github.com/JuliaArrays/FixedSizeArrays.jl/actions/workflows/UnitTests.yml/badge.svg?branch=main)](https://github.com/JuliaArrays/FixedSizeArrays.jl/actions/workflows/UnitTests.yml?query=branch%3Amain)
[![Coverage](https://codecov.io/gh/JuliaArrays/FixedSizeArrays.jl/branch/main/graph/badge.svg)](https://codecov.io/gh/JuliaArrays/FixedSizeArrays.jl)
[![PkgEval](https://JuliaCI.github.io/NanosoldierReports/pkgeval_badges/F/FixedSizeArrays.svg)](https://JuliaCI.github.io/NanosoldierReports/pkgeval_badges/F/FixedSizeArrays.html)

[`FixedSizeArrays.jl`](https://github.com/JuliaArrays/FixedSizeArrays.jl) is a package for the [Julia programming language](https://julialang.org/) which implements mutable fixed-size arrays (`FixedSizeArray`), which means the length of the array can not change after construction and is amenable to be [constant-propagated](https://en.wikipedia.org/wiki/Constant_folding) at compile-time when possible.

FixedSizeArrays supports the [standard array interfaces](https://docs.julialang.org/en/v1/manual/interfaces/#man-interface-array), so things like broadcasting, matrix multiplication, other linear algebra operations, `similar`, `copyto!` or `map` should just work.

Find more in the manual:

* [documentation of latest stable version](https://juliaarrays.github.io/FixedSizeArrays.jl/stable)
* [documentation of in-development version](https://juliaarrays.github.io/FixedSizeArrays.jl/dev)
