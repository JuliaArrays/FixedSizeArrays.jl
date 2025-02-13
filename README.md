## `FixedSizeArrays.jl`

[![Aqua QA](https://raw.githubusercontent.com/JuliaTesting/Aqua.jl/master/badge.svg)](https://github.com/JuliaTesting/Aqua.jl)
[![Build Status](https://github.com/JuliaArrays/FixedSizeArrays.jl/actions/workflows/UnitTests.yml/badge.svg?branch=main)](https://github.com/JuliaArrays/FixedSizeArrays.jl/actions/workflows/UnitTests.yml?query=branch%3Amain)
[![Coverage](https://codecov.io/gh/JuliaArrays/FixedSizeArrays.jl/branch/main/graph/badge.svg)](https://codecov.io/gh/JuliaArrays/FixedSizeArrays.jl)
[![PkgEval](https://JuliaCI.github.io/NanosoldierReports/pkgeval_badges/F/FixedSizeArrays.svg)](https://JuliaCI.github.io/NanosoldierReports/pkgeval_badges/F/FixedSizeArrays.html)

`FixedSizeArrays.jl` is a proof-of-concept package for the [Julia programming language](https://julialang.org/) which implements mutable fixed-size arrays, which means the length of the array is constant and is amenable to be [constant-propagated](https://en.wikipedia.org/wiki/Constant_folding) at compile-time when possible.
This is an alternative implementation to [`MArray`](https://juliaarrays.github.io/StaticArrays.jl/stable/pages/api/#StaticArraysCore.MArray) from [`StaticArrays.jl`](https://github.com/JuliaArrays/StaticArrays.jl).

Main differences between `FixedSizeArray` and `MArray` are:

* `FixedSizeArray` is based on the `Memory` type introduced in Julia v1.11, `MArray` is backed by tuples;
* the size of the array is part of the type parameters of `MArray`, this isn't the case for `FixedSizeArray`, where the size is only a constant field of the data structure.

FixedSizeArrays supports the usual array interfaces, so things like broadcasting, matrix
multiplication, other linear algebra operations, `similar`, `copyto!` or `map` should just work.

Use the constructors to convert from other array types. Use `collect_as` to convert from
arbitrary iterators.

```julia-repl
julia> arr = [10 20; 30 14]
2×2 Matrix{Int64}:
 10  20
 30  14

julia> iter = (i for i ∈ 7:9 if i≠8);

julia> using FixedSizeArrays

julia> FixedSizeArray(arr)  # construct from an `AbstractArray` value
2×2 FixedSizeMatrix{Int64}:
 10  20
 30  14

julia> FixedSizeArray{Float64}(arr)  # construct from an `AbstractArray` value while converting element type
2×2 FixedSizeMatrix{Float64}:
 10.0  20.0
 30.0  14.0

julia> const ca = FixedSizeArrays.collect_as
collect_as (generic function with 1 method)

julia> ca(FixedSizeArray, iter)  # construct from an arbitrary iterator
2-element FixedSizeVector{Int64}:
 7
 9

julia> ca(FixedSizeArray{Float64}, iter)  # construct from an arbitrary iterator while converting element type
2-element FixedSizeVector{Float64}:
 7.0
 9.0
```

Note: `FixedSizeArray`s are not guaranteed to be stack-allocated, in fact they will more likely *not* be stack-allocated.
However, in some *extremely* simple cases the compiler may be able to completely elide their allocations:
```julia
julia> using FixedSizeArrays

julia> @noinline f(A::AbstractArray) = length(A)
f (generic function with 1 method)

julia> g() = f(FixedSizeVector{Float64}(undef, 3))
g (generic function with 1 method)

julia> h() = f(Vector{Float64}(undef, 3))
h (generic function with 1 method)

julia> code_llvm(g)
```
```llvm
; Function Signature: g()
;  @ REPL[3]:1 within `g`
define i64 @julia_g_511() #0 {
top:
  ret i64 3
}

```
```julia
julia> code_llvm(h)
```
```llvm
; Function Signature: h()
;  @ REPL[4]:1 within `h`
define i64 @julia_h_693() #0 {
top:
  %gcframe1 = alloca [3 x ptr], align 16
  call void @llvm.memset.p0.i64(ptr align 16 %gcframe1, i8 0, i64 24, i1 true)
  %pgcstack = call ptr inttoptr (i64 7452881148 to ptr)(i64 262) #10
  store i64 4, ptr %gcframe1, align 16
  %task.gcstack = load ptr, ptr %pgcstack, align 8
  %frame.prev = getelementptr inbounds ptr, ptr %gcframe1, i64 1
  store ptr %task.gcstack, ptr %frame.prev, align 8
  store ptr %gcframe1, ptr %pgcstack, align 8
; ┌ @ boot.jl:576 within `Array`
; │┌ @ boot.jl:514 within `GenericMemory`
    %"Memory{Float64}[]" = call ptr @jl_alloc_genericmemory(ptr nonnull @"+Core.GenericMemory#695.jit", i64 3)
; │└
; │ @ boot.jl:577 within `Array`
   %.data_ptr = getelementptr inbounds { i64, ptr }, ptr %"Memory{Float64}[]", i64 0, i32 1
   %0 = load ptr, ptr %.data_ptr, align 8
   %gc_slot_addr_0 = getelementptr inbounds ptr, ptr %gcframe1, i64 2
   store ptr %"Memory{Float64}[]", ptr %gc_slot_addr_0, align 16
   %ptls_field = getelementptr inbounds ptr, ptr %pgcstack, i64 2
   %ptls_load = load ptr, ptr %ptls_field, align 8
   %"new::Array" = call noalias nonnull align 8 dereferenceable(32) ptr @ijl_gc_pool_alloc_instrumented(ptr %ptls_load, i32 800, i32 32, i64 4645053728) #8
   %"new::Array.tag_addr" = getelementptr inbounds i64, ptr %"new::Array", i64 -1
   store atomic i64 4645053728, ptr %"new::Array.tag_addr" unordered, align 8
   %1 = getelementptr inbounds ptr, ptr %"new::Array", i64 1
   store ptr %0, ptr %"new::Array", align 8
   store ptr %"Memory{Float64}[]", ptr %1, align 8
   %"new::Array.size_ptr" = getelementptr inbounds i8, ptr %"new::Array", i64 16
   store i64 3, ptr %"new::Array.size_ptr", align 8
   store ptr %"new::Array", ptr %gc_slot_addr_0, align 16
; └
  %2 = call i64 @j_f_699(ptr nonnull %"new::Array")
  %frame.prev10 = load ptr, ptr %frame.prev, align 8
  store ptr %frame.prev10, ptr %pgcstack, align 8
  ret i64 %2
}
```

> [!WARNING]
> This package should currently be used only to experiment with the idea of `Memory`-backed fixed-size arrays, it's highly non-optimised, absolutely don't use it for production.

### What about the other package with same name?

FixedSizeArrays.jl is not related to a package with the same name by @SimonDanisch.
That, earlier, package was one of the StaticArrays.jl-like packages in the pre-v1
days of Julia: https://github.com/SimonDanisch/FixedSizeArrays.jl
