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
# `fill!` for `FixedSizeArray`s and their views, and so does not vectorize the loop.
const fc_fill_bounds_check_skipped = VERSION >= v"1.13"

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
        # Do not print the FileCheck diagnostic while the test is known to be broken.
        dump_input = fc_fill_bounds_check_skipped ? "never" : "fail"
        @testset "$(name)" for (name, A) in (("array", fc_fill_array), ("view", fc_fill_view))
            @test (@filecheck dump_input=dump_input begin
                @check "define"
                @check_not "{{boundserror|bounds_error}}"
                @check "store <{{[0-9]+}} x double>"
                code_llvm(fill!, (A, Float64); debuginfo=:none)
            end) skip=fc_fill_bounds_check_skipped
        end
    end
end
