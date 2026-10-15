# FileCheck-based tests: assert the shape of the LLVM IR generated for small functions using
# `FixedSizeArray`s.  Sources for most examples:
# https://github.com/JuliaArrays/FixedSizeArrays.jl/discussions/62

using FileCheck
using InteractiveUtils: code_llvm

# FileCheck regex matching any runtime heap-allocation entry point, across Julia versions.
const fc_alloc_call = "{{(i?jl_gc_(small|pool|big)_alloc|i?jl_alloc_(genericmemory|array))}}"
# `Int` is `i32` on 32-bit platforms.
const fc_int = "i$(Sys.WORD_SIZE)"
# FileCheck regex matching a call to throw a bounds error, across Julia versions.
const fc_bounds_error = "{{boundserror|bounds_error}}"

# Elimination of the heap allocation of the underlying `Memory` only works reliably from
# Julia 1.12 on, and requires bounds checking not to be forced on.
const fc_run_alloc_elimination_checks = VERSION >= v"1.12" && Base.JLOptions().check_bounds != 1

# Proving that `reshape` to an incompatible size always throws worked on Julia 1.12-DEV in
# 2024, regressed in the 1.12 and 1.13 releases, and works again on 1.14-DEV.
const fc_run_reshape_throw_checks = VERSION >= v"1.14-DEV"

# Bounds checks can be elided only when bounds checking is not forced on.
const fc_run_bounds_check_elision_checks = Base.JLOptions().check_bounds != 1
# https://github.com/JuliaArrays/FixedSizeArrays.jl/issues/200: since Julia 1.13 (which
# dropped `@inbounds` from the generic `fill!`) LLVM does not elide the bounds check in
# `fill!` for views of `FixedSizeArray`s, and so does not vectorize the loop.
# `FixedSizeArray`s themselves have a `fill!` method of their own, which does not have
# this problem.  (`v"1.13-"` to include the prereleases of v1.13.)
const fc_fill_view_bounds_check_skipped = VERSION >= v"1.13-"

# Before Julia 1.13 `copy` went through a loop with a bounds check and two allocations.
const fc_run_copy_checks = fc_run_bounds_check_elision_checks && VERSION >= v"1.13"

# https://github.com/JuliaArrays/FixedSizeArrays.jl/issues/201: the dimensions of an
# argument must be loaded once, in the check of its size, and not again in the loops that
# follow.  The LLVM names of the loads (`%"[2]_ptr`) match from Julia 1.12 on.
const fc_run_size_reload_checks = VERSION >= v"1.12"

# Constant folding of `length` through a `@noinline` call.
@noinline fc_len(A::AbstractArray) = length(A)
fc_len_fsv() = fc_len(FixedSizeVector{Float64}(undef, 3))

# Constant folding of `length` after a side-effecting call.
function fc_rand_len()
    v = FixedSizeVector{Float64}(undef, 8)
    Random.rand!(v)
    return length(v)
end

# Reduction over a small random vector.
function fc_rand_sum()
    v = FixedSizeVector{Float64}(undef, 15)
    Random.rand!(v)
    return sum(v)
end

# Reduction over a larger random vector.
function fc_rand_foldl()
    v = FixedSizeVector{Float64}(undef, 250)
    Random.rand!(v)
    return foldl(+, v)
end

# Escape analysis removes the heap allocation even without `@inbounds`.
function fc_fill_reduce(typ::Type{<:AbstractArray}, n, d)
    @inline let
        a = typ{Float32}(undef, 15)
        for i ∈ eachindex(a)
            a[i] = n + i
        end
        a[d] = 0
        reduce(+, a)
    end
end
fc_fill_reduce_fsv(n, d) = @inline fc_fill_reduce(FixedSizeVectorDefault, n, d)

# Construction of an empty vector.
fc_empty() = FixedSizeVector{Float64}(undef, 0)

# Non-escaping vector of dynamic length.
fc_dynamic_len(n) = length(FixedSizeVector{Float64}(undef, n))

# `reshape` to a size incompatible with the (constant) length of the vector: the compiler
# can prove the call always throws.
@noinline fc_reshape(v::AbstractVector{Float64}) = reshape(v, 2, 2)
fc_reshape_3() = fc_reshape(FixedSizeVector{Float64}(undef, 3))
fc_reshape_4() = fc_reshape(FixedSizeVector{Float64}(undef, 4))

# `fill!` on an array and on a view along its last dimension, see issue #200.
const fc_fill_array = FixedSizeArrayDefault{Float64, 5}
const fc_fill_view = typeof(view(fc_fill_array(undef, 1, 1, 1, 1, 1), :, :, :, :, 1))

# Scalar indexing.
fc_getindex(A, i) = A[i]
fc_setindex!(A, i, x) = (A[i] = x; nothing)

# Cartesian indexing of an array with static dimensions.
const fc_static_dims_array = FixedSizeArrayDims{Float64, 3, Tuple{StaticInt{5}, StaticInt{4}, Int}, FixedSizeArrays.default_underlying_storage_type{Float64}}
fc_getindex(A, i, j, k) = A[i, j, k]
fc_getindex_inbounds(A, i, j, k) = @inbounds A[i, j, k]
fc_setindex_inbounds!(A, i, j, k, x) = (@inbounds A[i, j, k] = x; nothing)

# Loop over `eachindex`.
function fc_sum_loop(A)
    s = zero(eltype(A))
    for i in eachindex(A)
        s += A[i]
    end
    return s
end

# In-place broadcasting.
fc_broadcast!(A) = (A .= A .+ A; nothing)

# Constructor with constant dimensions.
fc_constant_dims() = FixedSizeMatrixDefault{Float64}(undef, 3, 4)

# `size` and `length` of an argument.
fc_size(A) = size(A)
fc_length(A) = length(A)

# `size` in a static dimension of an argument.
fc_static_size(A) = size(A, 2) * 3

# `copy`.
fc_copy(A) = copy(A)

# Reduced version of the kernel in issue #201: a check of the size of the argument
# followed by loops, one of which is conditionally executed, writing to `du`.  All three
# ingredients are needed to reproduce the issue: with a plain loop nest, or without the
# stores, or with two dimensions only, the IR is the same with lazy and eager `size`.
@noinline fc_201_throw_dm(ax) = throw(DimensionMismatch("unexpected axes $(ax)"))
function fc_201_kernel!(du, nelem)
    axes(du) == (Base.OneTo(5), Base.OneTo(4), Base.OneTo(nelem)) || fc_201_throw_dm(axes(du))
    for e in 1:nelem, i in 1:4
        for ii in (i + 1):4, v in 1:5
            @inbounds du[v, ii, e] += 1.0
        end
    end
    return nothing
end

# Broadcasting, whose aliasing check must use `Base.dataids` instead of `objectid`, see
# issues #19 and #63, and PR #64.
fc_broadcast_scale!(B, a, A) = (B .= a .* A; nothing)
fc_broadcast_copy!(B, A) = (B .= A; nothing)

# Out-of-bounds Cartesian indexing must not capture the array in the error, otherwise the
# heap allocation of the `Memory` cannot be eliminated, see issue #175 and PR #176.
function fc_fill_reduce_2d(n, i, j)
    a = FixedSizeMatrixDefault{Float64}(undef, 3, 4)
    fill!(a, n)
    a[i, j] = 0
    return sum(a)
end

# Stores in a loop without `@inbounds`, see issue #70 and PR #68.
function fc_store_loop!(v)
    for i in eachindex(v)
        v[i] = i
    end
    return nothing
end

# `@simd ivdep` loop over `eachindex` of two arrays, see issue #152.
const fc_coeffs = (0.1, 0.2, 0.3, 0.4, 0.5)
function fc_simd_ivdep!(out, x, c)
    @simd ivdep for i in eachindex(out, x)
        out[i] = muladd(c, evalpoly(x[i], fc_coeffs), out[i])
    end
    return nothing
end

# Iteration with `for x in v`, see PR #91.
function fc_iterate_sum(v)
    s = zero(eltype(v))
    for x in v
        s += x
    end
    return s
end

# Cartesian indexing.
fc_getindex2(A, i, j) = A[i, j]

# Splatting and destructuring, which need constant propagation through `iterate`, see
# PR #158.
function fc_splat()
    v = FixedSizeVector{Int}(undef, 3)
    v[1] = 1; v[2] = 2; v[3] = 3
    return tuple(v...)
end
function fc_destructure(x)
    v = FixedSizeVectorDefault{Float64}(undef, 3)
    v[1] = x; v[2] = 2x; v[3] = 3x
    a, b, c = v
    return a + b + c
end

# Destructuring with a slurp, which uses the `Base.rest` method of PR #166.
fc_rest(v) = ((_, r...) = v; r)

# Constant folding through `collect_as`, see issues #71 and #160, and PR #172.
fc_collect_as_len() = length(collect_as(FixedSizeVector{Int}, (1, 2, 3)))
fc_collect_as_sum() = sum(collect_as(FixedSizeVectorDefault{Int}, (1, 2, 3)))

# Constant folding of `length` with non-`Int` dimensions, see PR #47.
fc_len_int32() = fc_len(FixedSizeVector{Float64}(undef, Int32(3)))
fc_len_mixed() = fc_len(FixedSizeMatrix{Float64}(undef, Int32(3), 0x2))

# `isassigned`, see issue #14 and PR #15.
fc_isassigned(A, i) = isassigned(A, i)

# `A[:]` (PRs #130 and #161) and the copying constructor (issue #43).
fc_colon(A) = A[:]
fc_copy_constructor(A) = FixedSizeArray(A)

# `copy!` into a `FixedSizeVector`, which used to fall back to `resize!`, see issue #194 and
# PR #195.
fc_copy!(dst, src) = copy!(dst, src)

# The IR of a build instrumented for code coverage is not representative.
build_is_production_build &&
@testset "FileCheck: LLVM IR properties" begin
    @testset "`length` is constant-folded through a `@noinline` call" begin
        # The whole function body is a single `ret` of a literal.
        @test @filecheck begin
            @check "define"
            @check_next "top:"
            @check_next "ret $fc_int 3"
            @check_next "}"
            code_llvm(fc_len_fsv, (); debuginfo=:none)
        end
    end

    @testset "`length` is constant-folded after a side-effecting call" begin
        @test @filecheck begin
            @check "define"
            @check "ret $fc_int 8"
            code_llvm(fc_rand_len, (); debuginfo=:none)
        end
    end

    @testset "constructing an empty vector emits no allocation call" begin
        @test @filecheck begin
            @check "define"
            @check_not fc_alloc_call
            @check "ret void"
            code_llvm(fc_empty, (); debuginfo=:none)
        end
    end

    fc_run_alloc_elimination_checks &&
    @testset "heap allocation elimination" begin
        @testset "`length` after a side-effecting call" begin
            # The storage lives on the stack: no GC frame, no heap allocation call.
            @test @filecheck begin
                @check "define"
                @check_not "gcframe"
                @check_not fc_alloc_call
                @check "= alloca"
                @check_not fc_alloc_call
                @check "ret $fc_int 8"
                code_llvm(fc_rand_len, (); debuginfo=:none)
            end
        end

        @testset "sum of a small random vector: unrolled reduction" begin
            @test @filecheck begin
                @check "define"
                @check_not "gcframe"
                @check_not fc_alloc_call
                @check "= alloca"
                @check_not fc_alloc_call
                @check "fadd double"
                @check "ret double"
                code_llvm(fc_rand_sum, (); debuginfo=:none)
            end
        end

        @testset "`foldl` over a larger random vector" begin
            @test @filecheck begin
                @check "define"
                @check_not "gcframe"
                @check_not fc_alloc_call
                @check "= alloca"
                @check_not fc_alloc_call
                @check "ret double"
                code_llvm(fc_rand_foldl, (); debuginfo=:none)
            end
        end

        @testset "escape analysis without `@inbounds`" begin
            @test @filecheck begin
                @check "define"
                @check_not "gcframe"
                @check_not fc_alloc_call
                @check "= alloca"
                @check_not fc_alloc_call
                @check "ret float"
                code_llvm(fc_fill_reduce_fsv, (Int, Int); debuginfo=:none)
            end
        end

        @testset "non-escaping vector of dynamic length" begin
            # The allocation is removed and the length argument is returned directly.
            @test @filecheck begin
                @check "define"
                @check_not "gcframe"
                @check_not fc_alloc_call
                @check "ret $fc_int %\"n::Int$(Sys.WORD_SIZE)\""
                code_llvm(fc_dynamic_len, (Int,); debuginfo=:none)
            end
        end

        @testset "out-of-bounds Cartesian indexing (issue #175)" begin
            # A bounds error capturing the array would force the allocation of its `Memory`.
            @test @filecheck begin
                @check "define"
                @check_not "gcframe"
                @check_not fc_alloc_call
                @check "ret double"
                code_llvm(fc_fill_reduce_2d, (Float64, Int, Int); debuginfo=:none)
            end
        end

        @testset "splatting and destructuring are constant-folded (PR #158)" begin
            @test @filecheck begin
                @check "define"
                @check_not "gcframe"
                @check_not fc_alloc_call
                @check "store $fc_int 1"
                @check "store $fc_int 2"
                @check "store $fc_int 3"
                @check_not fc_alloc_call
                @check "ret void"
                code_llvm(fc_splat, (); debuginfo=:none)
            end
            @test @filecheck begin
                @check "define"
                @check_not "{{gcframe| call }}"
                @check "ret double"
                code_llvm(fc_destructure, (Float64,); debuginfo=:none)
            end
        end

        @testset "`collect_as` from a tuple is constant-folded (PR #172)" begin
            @testset "$(f)" for (f, result) in ((fc_collect_as_len, 3), (fc_collect_as_sum, 6))
                @test @filecheck begin
                    @check "define"
                    @check_not "{{ call |gcframe| br }}"
                    @check "ret $fc_int $result"
                    @check_next "}"
                    code_llvm(f, (); debuginfo=:none)
                end
            end
        end
    end

    fc_run_reshape_throw_checks &&
    @testset "`reshape` to an incompatible size is known to always throw" begin
        # The function ends with a trap instead of returning.
        @test @filecheck begin
            @check "define"
            @check "llvm.trap"
            @check_next "unreachable"
            @check_next "}"
            code_llvm(fc_reshape_3, (); debuginfo=:none)
        end
        # With a compatible size the function returns normally.
        @test @filecheck begin
            @check "define"
            @check_not "llvm.trap"
            @check "ret void"
            code_llvm(fc_reshape_4, (); debuginfo=:none)
        end
    end

    @testset "constructor with constant dimensions has no error branches" begin
        @test @filecheck begin
            @check "define"
            @check_not "{{throw|argument_error|unreachable}}"
            @check "ret void"
            code_llvm(fc_constant_dims, (); debuginfo=:none)
        end
    end

    @testset "`size` and `length` of an argument are straight-line loads" begin
        # A single basic block: no branches, no calls (so no `Memory` null check or length load).
        @test @filecheck begin
            @check "define"
            @check_next "top:"
            @check_not "{{ br | call }}"
            @check "ret void"
            code_llvm(fc_size, (FixedSizeMatrixDefault{Float64},); debuginfo=:none)
        end
        @test @filecheck begin
            @check "define"
            @check_next "top:"
            @check_not "{{ br | call }}"
            @check "ret $fc_int"
            code_llvm(fc_length, (FixedSizeMatrixDefault{Float64},); debuginfo=:none)
        end
    end

    @testset "`size` in a static dimension is constant-folded" begin
        # The whole function body is a single `ret` of a literal: the size is not loaded.
        @test @filecheck begin
            @check "define"
            @check_next "top:"
            @check_next "ret $fc_int 12"
            @check_next "}"
            code_llvm(fc_static_size, (fc_static_dims_array,); debuginfo=:none)
        end
    end

    fc_run_size_reload_checks &&
    @testset "dimensions are not reloaded after a check of the size (issue #201)" begin
        # With `size` loading each dimension lazily, the second dimension of `du` was
        # reloaded in the loop and used as a dynamic stride.  The only load of the second
        # dimension must be the one in the size check, and the strides must all be
        # constants: no multiplication by a non-constant value.
        @test @filecheck begin
            @check "define"
            @check_not "mul $fc_int %{{[^,]+}}, %"
            @check "load $fc_int, ptr %\"[2]_ptr"
            @check_not "load $fc_int, ptr %\"[2]_ptr"
            @check_not "mul $fc_int %{{[^,]+}}, %"
            @check "ret void"
            code_llvm(fc_201_kernel!, (FixedSizeArrayDefault{Float64, 3}, Int); debuginfo=:none)
        end
    end

    # On Julia 1.10 the call throwing the bounds error goes through an anonymous pointer, so
    # it cannot be matched by name.
    fc_run_bounds_check_elision_checks && VERSION >= v"1.11" &&
    @testset "scalar indexing performs exactly one bounds check" begin
        # The `@inbounds` access to the parent must not add a second check.
        @testset "`getindex`" begin
            @test @filecheck begin
                @check "define"
                @check_count 1 fc_bounds_error
                @check_not fc_bounds_error
                code_llvm(fc_getindex, (FixedSizeVectorDefault{Float64}, Int); debuginfo=:none)
            end
        end
        @testset "`setindex!`" begin
            @test @filecheck begin
                @check "define"
                @check_count 1 fc_bounds_error
                @check_not fc_bounds_error
                code_llvm(fc_setindex!, (FixedSizeVectorDefault{Float64}, Int, Float64); debuginfo=:none)
            end
        end
    end

    fc_run_bounds_check_elision_checks && VERSION >= v"1.11" &&
    @testset "Cartesian indexing of an array with static dimensions" begin
        @testset "one bounds check" begin
            @test @filecheck begin
                @check "define"
                @check_count 1 fc_bounds_error
                @check_not fc_bounds_error
                code_llvm(fc_getindex, (fc_static_dims_array, Int, Int, Int); debuginfo=:none)
            end
        end
        # The static sizes are constants (e.g. the multiplication by the first size), and the
        # size in the last dimension is not loaded.
        @testset "`@inbounds`: no bounds check, static sizes" begin
            @test @filecheck begin
                @check "define"
                @check_not "{{bounds_error|boundserror|load i64}}"
                @check "mul i64 {{.*}}, 5"
                @check_not "{{bounds_error|boundserror|load i64}}"
                @check "load double"
                @check_not "{{bounds_error|boundserror|load i64}}"
                @check "ret double"
                code_llvm(fc_getindex_inbounds, (fc_static_dims_array, Int, Int, Int); debuginfo=:none)
            end
            @test @filecheck begin
                @check "define"
                @check_not fc_bounds_error
                @check "store double"
                @check_not fc_bounds_error
                @check "ret void"
                code_llvm(fc_setindex_inbounds!, (fc_static_dims_array, Int, Int, Int, Float64); debuginfo=:none)
            end
        end
    end

    fc_run_bounds_check_elision_checks &&
    @testset "loop over `eachindex` has no bounds check" begin
        @test @filecheck begin
            @check "define"
            @check_not fc_bounds_error
            @check "ret double"
            code_llvm(fc_sum_loop, (FixedSizeVectorDefault{Float64},); debuginfo=:none)
        end
    end

    fc_run_bounds_check_elision_checks &&
    @testset "in-place broadcasting: no allocation, no bounds check, vectorized" begin
        @test @filecheck begin
            @check "define"
            @check_not fc_alloc_call
            @check_not fc_bounds_error
            @check "fadd <{{[0-9]+}} x double>"
            @check_not fc_alloc_call
            @check_not fc_bounds_error
            @check "ret void"
            code_llvm(fc_broadcast!, (FixedSizeVectorDefault{Float64},); debuginfo=:none)
        end
        @test @filecheck begin
            @check "define"
            @check_not fc_alloc_call
            @check_not fc_bounds_error
            @check "fadd <{{[0-9]+}} x double>"
            @check_not fc_alloc_call
            @check_not fc_bounds_error
            @check "ret void"
            code_llvm(fc_broadcast!, (fc_static_dims_array,); debuginfo=:none)
        end
    end

    fc_run_copy_checks &&
    @testset "`copy` is one allocation and a bulk memory copy" begin
        @test @filecheck begin
            @check "define"
            @check_not "{{vector.body|boundserror|bounds_error}}"
            @check "jl_alloc_genericmemory"
            @check_not "{{vector.body|boundserror|bounds_error}}"
            @check "memmove"
            @check_not "{{vector.body|boundserror|bounds_error}}"
            @check "ret void"
            code_llvm(fc_copy, (FixedSizeVectorDefault{Float64},); debuginfo=:none)
        end
    end

    fc_run_bounds_check_elision_checks &&
    @testset "`fill!` elides the bounds check and vectorizes (issue #200)" begin
        @testset "$(name)" for (name, A, skipped) in (
            ("array", fc_fill_array, false),
            ("static dimensions", fc_static_dims_array, false),
            ("view", fc_fill_view, fc_fill_view_bounds_check_skipped),
        )
            # Do not print the FileCheck diagnostic while the test is known to be broken.
            dump_input = skipped ? "never" : "fail"
            @test (@filecheck dump_input=dump_input begin
                @check "define"
                @check_not "{{boundserror|bounds_error}}"
                @check "store <{{[0-9]+}} x double>"
                code_llvm(fill!, (A, Float64); debuginfo=:none)
            end) skip=skipped
        end
    end

    @testset "broadcasting does not call `objectid` (issue #63)" begin
        @test @filecheck begin
            @check "define"
            @check_not "object_id"
            code_llvm(fc_broadcast_scale!, (FixedSizeVectorDefault{Float64}, Float64, FixedSizeVectorDefault{Float64}); debuginfo=:none)
        end
        @test @filecheck begin
            @check "define"
            @check_not "object_id"
            code_llvm(fc_broadcast_copy!, (FixedSizeVectorDefault{Float64}, FixedSizeVectorDefault{Float64}); debuginfo=:none)
        end
    end

    fc_run_bounds_check_elision_checks && VERSION >= v"1.13" &&
    @testset "loops without `@inbounds` are vectorized" begin
        @testset "stores (issue #70)" begin
            @test @filecheck begin
                @check "define"
                @check_not fc_bounds_error
                @check "store <{{[0-9]+}} x double>"
                @check_not fc_bounds_error
                code_llvm(fc_store_loop!, (FixedSizeVectorDefault{Float64},); debuginfo=:none)
            end
        end
        @testset "`@simd ivdep` over `eachindex` of two arrays (issue #152)" begin
            @test @filecheck begin
                @check "define"
                @check_not fc_bounds_error
                @check "vector.body"
                @check "<{{[0-9]+}} x double>"
                @check_not fc_bounds_error
                code_llvm(fc_simd_ivdep!, (FixedSizeVectorDefault{Float64}, FixedSizeVectorDefault{Float64}, Float64); debuginfo=:none)
            end
        end
    end

    fc_run_bounds_check_elision_checks &&
    @testset "loop over `eachindex` of an N-d array: no bounds check, no division (PR #6)" begin
        # Linear indexing: no conversion of the index to Cartesian indices.
        @testset "N = $N" for N in 1:3
            @test @filecheck begin
                @check "define"
                @check_not "{{boundserror|bounds_error|sdiv|udiv|srem|urem}}"
                @check "ret double"
                code_llvm(fc_sum_loop, (FixedSizeArrayDefault{Float64, N},); debuginfo=:none)
            end
        end
    end

    fc_run_bounds_check_elision_checks &&
    @testset "`for x in v` has no bounds check and is vectorized (PR #91)" begin
        @test @filecheck begin
            @check "define"
            @check_not fc_bounds_error
            @check "add <{{[0-9]+}} x i8>"
            @check_not fc_bounds_error
            @check "ret i8"
            code_llvm(fc_iterate_sum, (FixedSizeVectorDefault{Int8},); debuginfo=:none)
        end
    end

    # On Julia 1.10 and 1.11 the bounds error is boxed on the heap.
    fc_run_bounds_check_elision_checks && VERSION >= v"1.12" &&
    @testset "Cartesian indexing performs exactly one bounds check and no allocation" begin
        @test @filecheck begin
            @check "define"
            @check_not fc_alloc_call
            @check_count 1 fc_bounds_error
            @check_not fc_bounds_error
            @check_not fc_alloc_call
            code_llvm(fc_getindex2, (FixedSizeMatrixDefault{Float64}, Int, Int); debuginfo=:none)
        end
    end

    @testset "splatting does not go through dynamic dispatch (PR #158)" begin
        @test @filecheck begin
            @check "define"
            @check_not "{{apply_iterate|jl_apply_generic|jl_invoke}}"
            @check "ret"
            code_llvm(fc_splat, (); debuginfo=:none)
        end
    end

    @testset "destructuring with a slurp uses the specialized `Base.rest` (PR #166)" begin
        @test @filecheck begin
            @check "define"
            @check_not "{{_foldl_impl|sizehint}}"
            code_llvm(fc_rest, (FixedSizeVectorDefault{Float64},); debuginfo=:none)
        end
    end

    @testset "`length` is constant-folded with non-`Int` dimensions (PR #47)" begin
        @test @filecheck begin
            @check "define"
            @check_next "top:"
            @check_next "ret $fc_int 3"
            @check_next "}"
            code_llvm(fc_len_int32, (); debuginfo=:none)
        end
        # Julia 1.10 does not fold the product of mixed integer types.
        VERSION >= v"1.11" && @test @filecheck begin
            @check "define"
            @check_next "top:"
            @check_next "ret $fc_int 6"
            @check_next "}"
            code_llvm(fc_len_mixed, (); debuginfo=:none)
        end
    end

    fc_run_bounds_check_elision_checks &&
    @testset "`isassigned` is a bounds check and a null check (PR #15)" begin
        # No runtime calls.  Note that Julia 1.10 calls the `llvm.smax` intrinsic.
        @test @filecheck begin
            @check "define"
            @check_not "{{gcframe|jl_|j_}}"
            @check "icmp ne"
            @check_not "{{gcframe|jl_|j_}}"
            @check "ret i8"
            code_llvm(fc_isassigned, (FixedSizeVectorDefault{String}, Int); debuginfo=:none)
        end
    end

    fc_run_copy_checks &&
    @testset "`A[:]` and the copying constructor are one allocation and a bulk memory copy" begin
        @testset "$(f) $(T)" for (f, T) in (
                (fc_colon, FixedSizeVectorDefault{Float64}),
                (fc_colon, FixedSizeMatrixDefault{Float64}),
                (fc_copy_constructor, FixedSizeVectorDefault{Float64}),
            )
            @test @filecheck begin
                @check "define"
                @check_not "{{vector.body|boundserror|bounds_error}}"
                @check "jl_alloc_genericmemory"
                @check_not "{{vector.body|boundserror|bounds_error}}"
                @check "memmove"
                @check_not "{{vector.body|boundserror|bounds_error}}"
                @check "ret void"
                code_llvm(f, (T,); debuginfo=:none)
            end
        end
    end

    # On Julia 1.10 the storage is a `Vector` and the copy happens in a non-inlined call.
    fc_run_bounds_check_elision_checks && VERSION >= v"1.11" &&
    @testset "`copyto!` between `FixedSizeArray`s is a bulk memory copy" begin
        # Inspect `copyto!` directly: in a wrapper the call to it may not be inlined.
        @test @filecheck begin
            @check "define"
            @check_not "vector.body"
            @check "memmove"
            @check_not "vector.body"
            @check "ret"
            code_llvm(copyto!, (FixedSizeVectorDefault{Float64}, FixedSizeVectorDefault{Float64}); debuginfo=:none)
        end
    end

    @testset "`copy!` into a `FixedSizeVector` does not fall back to `resize!` (issue #194)" begin
        # On Julia 1.12 `copy!` is not inlined, so this does not catch the issue there.
        @test @filecheck begin
            @check "define"
            @check_not "{{apply_generic|methoderror}}"
            code_llvm(fc_copy!, (FixedSizeVectorDefault{Float64}, Vector{Float64}); debuginfo=:none)
        end
    end
end
