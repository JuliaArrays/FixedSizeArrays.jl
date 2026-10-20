using Test
using LinearAlgebra, Adapt, MPI
using FixedSizeArrays
using FixedSizeArrays: checked_dims, wrap, FixedSizeArrayDims
using Collects: collect_as
using OffsetArrays: OffsetArray
using Random: Random
using StaticArrayInterface: StaticArrayInterface
using Static: static, StaticInt
import Aqua

# The memory backends to test with.
const storage_types = (FixedSizeArrays.optional_memory..., Vector)

# Check if the compilation options allow maximum performance.
const build_is_production_build_env_name = "BUILD_IS_PRODUCTION_BUILD"
const build_is_production_build = let v = get(ENV, build_is_production_build_env_name, "true")
    if v ∉ ("false", "true")
        error("unknown value for environment variable $build_is_production_build_env_name: $v")
    end
    if v == "true"
        true
    else
        false
    end
end::Bool

# helpers for testing for allocation or suboptimal inference

allocated(f::F,      args::Tuple) where {F} = @allocated f(args...)
allocated(::Type{T}, args::Tuple) where {T} = @allocated T(args...)
test_noalloc(f::F,      args::Tuple) where {F} = build_is_production_build && @test iszero(allocated(f, args))
test_noalloc(::Type{T}, args::Tuple) where {T} = build_is_production_build && @test iszero(allocated(T, args))
function test_inferred(f::F,      ::Type{R}, args::Tuple) where {F,R}
    @test isconcretetype(R)  # meta: test-strictness test
    @test (@inferred f(args...)) isa R
end
function test_inferred(::Type{T}, ::Type{R}, args::Tuple) where {T,R}
    @test isconcretetype(R)  # meta: test-strictness test
    @test (@inferred T(args...)) isa R
end
test_inferred(::Type{T}, args::Tuple) where {T} = test_inferred(T, T, args)
function test_inferred_noalloc(f::F,      ::Type{R}, args::Tuple) where {F,R}
    test_noalloc(f, args)
    test_inferred(f, R, args)
end
function test_inferred_noalloc(::Type{T}, ::Type{R}, args::Tuple) where {T,R}
    test_noalloc(T, args)
    test_inferred(T, R, args)
end
function test_inferred_noalloc(::Type{T},            args::Tuple) where {T  }
    test_noalloc(T, args)
    test_inferred(T, args)
end

# helpers for constructing the type constructors

function fsa(vec_type::Type{<:DenseVector})
    FixedSizeArray{T,N,vec_type{T}} where {T,N}
end
function fsm(vec_type::Type{<:DenseVector})
    fsa(vec_type){T,2} where {T}
end
function fsv(vec_type::Type{<:DenseVector})
    fsa(vec_type){T,1} where {T}
end

function test_we_do_not_own_the_call(func, arg_types)
    function f(method)
        FixedSizeArrays != parentmodule(method)
    end
    @test all(f, methods(func, arg_types))
end


# helpers for the MPI tests

# `isbits` struct with padding, for which MPI.jl creates a derived datatype
struct MPITestStruct
    a::Float64
    b::Int32
end

mpi_test_value(::Type{T}, i::Integer) where {T<:Real} = T(i)
mpi_test_value(::Type{T}, i::Integer) where {T<:Complex} = T(i, -i)
mpi_test_value(::Type{MPITestStruct}, i::Integer) = MPITestStruct(i, -i)

# `DenseVector` which is neither a `Vector` nor a `Memory`, for testing that conversion to
# `MPI.MPIPtr` is delegated to the parent of the `FixedSizeArray` (as for GPU arrays)
struct MPITestVector{T} <: DenseVector{T}
    v::Vector{T}
end
Base.size(x::MPITestVector) = size(x.v)
Base.getindex(x::MPITestVector, i::Int) = x.v[i]
Base.setindex!(x::MPITestVector, y, i::Int) = (x.v[i] = y; x)
Base.cconvert(::Type{MPI.MPIPtr}, x::MPITestVector) = Base.cconvert(MPI.MPIPtr, x.v)

@testset verbose=true "FixedSizeArrays" begin
    @testset "meta" begin
        @test isempty(detect_ambiguities(Main))  # test for ambiguities in this file
    end

    @testset "Aqua.jl" begin
        Aqua.test_all(FixedSizeArrays)
    end

    @testset "type piracy" begin
        @testset "issue #77: type piracy of `similar`" begin
            test_we_do_not_own_the_call(similar, Tuple{Broadcast.Broadcasted{Broadcast.ArrayStyle{Union{}}}, Type{Int}})
        end
    end

    @testset "type aliases" begin
        @test FixedSizeArrayDefault <: FixedSizeArray
        @test FixedSizeVectorDefault <: FixedSizeVector
        @test FixedSizeMatrixDefault <: FixedSizeMatrix
        for et ∈ (Float32, Int, String)
            @test FixedSizeArrayDefault{et} <: FixedSizeArray{et}
            @test FixedSizeVectorDefault{et} <: FixedSizeVector{et}
            @test FixedSizeMatrixDefault{et} <: FixedSizeMatrix{et}
            @test FixedSizeArrayDefault{et, 0} <: FixedSizeArray{et, 0}
            for t ∈ (
                FixedSizeVectorDefault{et}, FixedSizeMatrixDefault{et}, FixedSizeArrayDefault{et, 0},
            )
                @test isconcretetype(t)
            end
        end
        @test FixedSizeArrays.FixedSizeArrayAllowedConstructorType <: Type{<:FixedSizeArray}
    end

    @testset "safe computation of length from dimensions size" begin
        @test isone(checked_dims(()))
        max_size = Sys.WORD_SIZE == 32 ? 12 : 20
        for n ∈ 0:30
            t = Tuple(1:n)
            if n > max_size
                @test_throws ArgumentError checked_dims(t)
            else
                @test factorial(n) == prod(t) == checked_dims(t)
                test_inferred_noalloc(checked_dims, Int, (t,))
            end
        end
    end

    struct Iter{E,N,I<:Integer}
        size::NTuple{N,I}
        length::I
        val::E
    end
    function Base.iterate(i::Iter)
        l = i.length
        iterate(i, max(zero(l), l))
    end
    function Base.iterate(i::Iter, state::Int)
        if iszero(state)
            nothing
        else
            (i.val, state - 1)
        end
    end
    Base.IteratorSize(::Type{<:Iter{<:Any,N}}) where {N} = Base.HasShape{N}()
    Base.length(i::Iter) = i.length
    Base.size(i::Iter) = i.size
    Base.eltype(::Type{<:Iter{E}}) where {E} = E

    @testset "default underlying storage type" begin
        default = FixedSizeArrays.default_underlying_storage_type
        @test default === (@isdefined(Memory) ? Memory : Vector)
        return_type = FixedSizeVector{Int,default{Int}}
        @test return_type === FixedSizeVectorDefault{Int}
        test_inferred(FixedSizeArray{Int}, return_type, (undef, 3))
        test_inferred(FixedSizeArray{Int}, return_type, (undef, (3,)))
        test_inferred(FixedSizeVector{Int}, return_type, (undef, 3))
        test_inferred(FixedSizeVector{Int}, return_type, (undef, (3,)))
        iter = Iterators.filter(iseven, 3:7)
        @test collect_as(FixedSizeArray, iter) isa return_type
        test_inferred(collect_as, return_type, (FixedSizeArray{Int}, iter))
        test_inferred(collect_as, return_type, (FixedSizeVector{Int}, iter))
        arr = ([1, 2, 3],)
        test_inferred(FixedSizeArray, return_type, arr)
        test_inferred(FixedSizeVector, return_type, arr)
        test_inferred(FixedSizeArray{Int}, return_type, arr)
        test_inferred(FixedSizeVector{Int}, return_type, arr)
    end

    @testset "`undef` construction without element type" begin
        @test_throws MethodError FixedSizeVector(undef, 1)
        @test_throws MethodError FixedSizeMatrix(undef, 1, 1)
    end

    # only run this testset when:
    #
    # * Julia is new enough for `Memory` to be available
    #
    # * bounds-checking is not forced on (so `--check-bounds=auto` or `--check-bounds=no`, but not `--check-bounds=yes`)
    (@isdefined Memory) && (Base.JLOptions().check_bounds != 1) &&
    @testset "examples of heap allocation being optimized out" begin
        # source for some/all of these: https://github.com/JuliaArrays/FixedSizeArrays.jl/discussions/62
        @testset "constant folding: 1" begin
            function f()
                v = FixedSizeVector{Float32}(undef, 9)
                @noinline length(v)
            end
            test_inferred_noalloc(f, Int, ())
        end
        @testset "constant folding: 2" begin
            function f()
                v = FixedSizeVector{Float32}(undef, 9)
                Random.rand!(v)
                length(v)
            end
            test_inferred_noalloc(f, Int, ())
        end
        @testset "constant folding: 3" begin
            function f()
                v = FixedSizeVector{Float32}(undef, 9)
                Random.rand!(v)
                sum(v)
            end
            test_inferred_noalloc(f, Float32, ())
        end
        @testset "constant folding: 4" begin
            function f()
                v = FixedSizeVector{Float32}(undef, 9)
                Random.rand!(v)
                foldl(+, v)
            end
            test_inferred_noalloc(f, Float32, ())
        end
        @testset "constant folding: 5" begin
            function f()
                v = FixedSizeVector{Float32}(undef, 9)
                Random.rand!(v)
                v[end]
            end
            test_inferred_noalloc(f, Float32, ())
        end
        @testset "allocation elimination without constant folding" begin
            function f(n, d)
                @inline let
                    a = FixedSizeVector{Float32}(undef, 9)
                    for i ∈ eachindex(a)
                        a[i] = n + i
                    end
                    a[d] = 0
                    reduce(+, a)
                end
            end
            test_inferred_noalloc(f, Float32, (123, 3))
        end
    end

    include("filecheck.jl")

    @testset verbose=true "test sets for multiple supported storage types" begin
        for storage_type ∈ storage_types
            FSV = fsv(storage_type)
            FSM = fsm(storage_type)
            FSA = fsa(storage_type)

            @testset "`propertynames`" begin
                for dim_count ∈ 0:4
                    for T ∈ (Int, Float32)
                        siz = ntuple(Returns(1), dim_count)
                        arr = FSA{T}(undef, siz)
                        test_inferred(propertynames, Tuple{}, (arr,))
                        test_inferred(propertynames, Tuple{}, (arr, false))
                    end
                end
            end

            @testset "Constructors" begin
                @testset "`undef` construction without element type" begin
                    @test_throws MethodError FSV(undef, 1)
                    @test_throws MethodError FSM(undef, 1, 1)
                end
                @testset "mismatched ndims" begin
                    @test_throws DimensionMismatch FSV{Int}(undef)
                    @test_throws DimensionMismatch FSV{Int}(undef, 1, 1)
                    @test_throws DimensionMismatch FSM{Int}(undef)
                    @test_throws DimensionMismatch FSM{Int}(undef, 1)
                end
                for dim_count ∈ 0:4
                    siz = ntuple(Returns(2), dim_count)
                    T = FSA{Float64}
                    R = FSA{Float64,dim_count}
                    for sz ∈ (siz, map(BigInt, siz))
                        for args ∈ ((undef, sz), (undef, sz...))
                            test_inferred(   R, args)
                            test_inferred(T, R, args)
                        end
                    end
                end
                for offset ∈ (-1, 0, 2, 3)
                    ax = offset:(offset + 1)
                    oa = OffsetArray([10, 20], ax)
                    @test_throws DimensionMismatch FSA{Int}(oa)
                    @test_throws DimensionMismatch FSV{Int}(oa)
                end
                @testset "taken from Julia's test/core.jl" begin
                    # inspired by:
                    # https://github.com/JuliaLang/julia/blob/83929ad883c97fa4376618c4559b74f6ada7a3ce/test/core.jl#L7211-L7228
                    b = prevpow(2, typemax(Int))
                    test_inferred(FSA{Int}, FSA{Int,3}, (undef, 0, b, b))
                    test_inferred(FSA{Int}, FSA{Int,3}, (undef, b, b, 0))
                    test_inferred(FSA{Int}, FSA{Int,5}, (undef, b, b, 0, b, b))
                    @test_throws ArgumentError FSA{Int}(undef, b, b)
                    @test_throws ArgumentError FSA{Int}(undef, 1, b, b)
                    @test_throws ArgumentError FSA{Int}(undef, 0, -10)
                    @test_throws ArgumentError FSA{Int}(undef, -10, 0)
                    @test_throws ArgumentError FSA{Int}(undef, -1, -1)
                    @test_throws ArgumentError FSA{Int}(undef, 0, -4, -4)
                    @test_throws ArgumentError FSA{Int}(undef, -4, 1, 0)
                    @test_throws ArgumentError FSA{Int}(undef, -4, -4, 1)
                end
                @test_throws ArgumentError FSA{Float64,1}(undef, -1)
                @test_throws ArgumentError FSA{Float64,1}(undef, (-1,))
                @test_throws ArgumentError FSA{Float64,2}(undef, -1, -1)
                @test_throws ArgumentError FSA{Float64,2}(undef, typemax(Int), 0)
                @test_throws ArgumentError FSA{Float64,2}(undef, typemax(Int), typemax(Int))
                @test_throws ArgumentError FSA{Float64,3}(undef, typemax(Int)-1, typemax(Int)-1, 2)
                @test_throws ArgumentError FSA{Float64,4}(undef, typemax(Int)-1, typemax(Int)-1, 2, 4)
                @testset "negative dimension size" begin
                    for n ∈ 2:4
                        funs = (
                            Returns(-1),
                            (i -> (i == 1) ?  1 : -1),
                            (i -> (i == 1) ? -1 :  1),
                            (i -> (i == n) ?  1 : -1),
                            (i -> (i == n) ? -1 :  1),
                        )
                        fun = f -> ntuple(f, n)
                        sizes = map(fun, funs)
                        for siz ∈ sizes
                            @test_throws ArgumentError FSA{Float64,n}(undef, siz)
                            @test_throws ArgumentError FSA{Float64,n}(undef, siz...)
                        end
                    end
                end
            end

            @testset "constructors taking `AbstractArray`" begin
                @test [1, 2, 3] == (@inferred FSV([1, 2, 3]))::FSV{Int}
                @test [1, 2, 3] == (@inferred FSA([1, 2, 3]))::FSV{Int}
                @test [1, 2, 3] == (@inferred FSV(Any[1, 2, 3]))::FSV{Any}
                @test [1, 2, 3] == (@inferred FSA(Any[1, 2, 3]))::FSV{Any}
                @test [1, 2, 3] == (@inferred FSV(Real[1, 2, 3]))::FSV{Real}
                @test [1, 2, 3] == (@inferred FSA(Real[1, 2, 3]))::FSV{Real}
                @test [1, 2, 3] == (@inferred FSV{Float32}([1, 2, 3]))::FSV{Float32}
                @test [1, 2, 3] == (@inferred FSA{Float32}([1, 2, 3]))::FSV{Float32}
                @test_throws DimensionMismatch FSV(rand(Bool, 2, 2))
                @test_throws DimensionMismatch FSV{Float32}(rand(Bool, 2, 2))
                @test let v = FSV([1, 2, 3])
                    v == (@inferred FSV(v))::typeof(v)
                end
                @test let v = FSV([1, 2, 3])
                    v !== FSV(v)
                end
            end

            @testset "setindex!" begin
                v = FSV{Float64}(undef, 2)
                v[1] = 1
                v[2] = 2
                @test v[1] == 1
                @test v[2] == 2
            end

            @testset "out-of-bounds access" begin
                v0 = FSV([])
                v1 = FSV([3])
                @test_throws BoundsErrorLight v0[1]
                @test_throws BoundsErrorLight v1[2]
                @test_throws BoundsErrorLight v0[1] = 7
                @test_throws BoundsErrorLight v1[2] = 7
                @test_throws ["LightBoundsError: ", "[2]"] v1[2]
                @test_throws BoundsErrorLight v0[]
                @test_throws BoundsErrorLight v0[1, 1]
                @test_throws BoundsErrorLight v0[1, 1, 1]
            end

            @testset "parent and parentindices" begin
                for size in ((4,), (3, 3), (2, 2, 2)), T in (UInt16, Int, Float32)
                    a = FixedSizeArray{T}(undef, size)
                    a .= 1
                    p = parent(a)
                    @test p === a.mem
                    @test parentindices(a) === (Base.OneTo(length(a)),)
                    # Make sure all the parent indices are valid
                    for idx in Base.CartesianIndex.(collect.(parentindices(a))...)
                        @test p[idx] == 1
                    end
                end
            end

            @testset "`FixedSizeArrays.wrap`" begin
                # `T` must be a static parameter here rather than a variable captured from the
                # `for T` loop below: a closure capturing `T` holds it as a `DataType`, so the
                # `wrap` call inside would be dynamically dispatched, and allocate (Julia v1.11).
                function test_wrap_noalloc(::Type{T}, p) where {T}
                    test_noalloc(p -> wrap(FixedSizeArray{T,3}, p, 2, 3, 4), (p,))
                end
                for T in (UInt16, Int, Float32)
                    p = storage_type{T}(undef, 24)
                    p .= 1:24
                    for size in ((24,), (4, 6), (2, 3, 4), (2, 3, 2, 2))
                        N = length(size)
                        R = FixedSizeArray{T,N,storage_type{T}}
                        for FSA_ in (FixedSizeArray, FixedSizeArray{T}, FixedSizeArray{T,N}, FSA, FSA{T}, FSA{T,N}, R)
                            a = wrap(FSA_, p, size)
                            @test a isa R
                            @test parent(a) === p
                            @test Base.size(a) === size
                            @test a == reshape(1:24, size)
                            @test wrap(FSA_, p, size...) === a
                            test_inferred(wrap, R, (FSA_, p, size))
                        end
                        # no copy: writes are visible in both directions
                        a = wrap(FixedSizeArray, p, size)
                        a[end] = 100
                        @test p[end] == 100
                        p[1] = 200
                        @test a[1] == 200
                        p .= 1:24
                    end
                    if storage_type !== Vector  # the `Vector` backend allocates on v1.10
                        test_noalloc(p -> wrap(FixedSizeArray, p, (4, 6)), (p,))
                        test_wrap_noalloc(T, p)
                    end
                    @test wrap(FixedSizeArray, p) isa FixedSizeArray{T,1,storage_type{T}}
                    @test Base.size(wrap(FixedSizeArray, p)) === (24,)
                    @test Base.size(wrap(FixedSizeArray, p, Int8(4), UInt(6))) === (4, 6)
                    # length mismatch
                    @test_throws DimensionMismatch wrap(FixedSizeArray, p, (5, 5))
                    @test_throws DimensionMismatch wrap(FixedSizeArray, p, 23)
                    @test_throws DimensionMismatch wrap(FixedSizeArray, p, ())
                    # invalid dimensions
                    @test_throws ArgumentError wrap(FixedSizeArray, p, (-4, -6))
                    # type parameters not matching `parent` or `size`
                    @test_throws ArgumentError wrap(FixedSizeArray{Complex{T}}, p, (4, 6))
                    @test_throws DimensionMismatch wrap(FixedSizeArray{T,3}, p, (4, 6))
                    @test_throws DimensionMismatch wrap(FixedSizeVector, p, (4, 6))
                    if @isdefined Memory
                        other_storage = storage_type === Vector ? Memory : Vector
                        @test_throws ArgumentError wrap(fsa(other_storage), p, (4, 6))
                    end
                end
                # 0-dimensional
                p = storage_type{Int}(undef, 1)
                p[1] = 7
                a = wrap(FixedSizeArray, p, ())
                @test a isa FixedSizeArray{Int,0,storage_type{Int}}
                @test a[] == 7
                # re-wrapping a `Vector` after `resize!`
                if storage_type === Vector
                    v = collect(1:6)
                    a = wrap(FixedSizeArray, v, (2, 3))
                    resize!(v, 8)
                    b = wrap(FixedSizeArray, v, (2, 4))
                    @test parent(b) === v
                    @test b == reshape(v, 2, 4)
                end
            end

            @testset "destructuring" begin
                vec = FSV(rand(3))
                mat = FSM(rand(3, 3))
                @testset "`Base.rest`" begin
                    @test let (v...,) = vec
                        v isa AbstractVector
                    end
                    @test let (v...,) = mat
                        v isa AbstractVector
                    end
                    @test let (v...,) = vec
                        v == vec
                    end
                    @test let (v...,) = mat
                        v == reshape(mat, :)
                    end
                end
            end

            @testset "FixedSizeVector" begin
                v = FSV{Float64}(undef, 3)
                @test length(v) == 3
                test_inferred_noalloc(((a, b) -> a .= b), FSV{Float64}, (v, 1:3))
                @test v == 1:3
                test_inferred(similar, FSV{Float64}, (v,))
                test_inferred(similar, FSV{Int}, (v, Int))
                test_inferred(copy, FSV{Float64}, (v,))
                test_inferred(zero, FSV{Float64}, (v,))
                test_inferred(similar, FSV{Int}, (FSV{Int}, (2,)))
                test_inferred(similar, FSV{Int}, (FSA{Int}, (2,)))
                example_abstract_vectors = (7:9, [7, 8, 9], OffsetArray([7, 8, 9], 1:3))
                for T in (FSA{Int}, FSV{Int}), arr in example_abstract_vectors
                    # Don't allow conversion:
                    # <https://github.com/JuliaArrays/FixedSizeArrays.jl/issues/120>.
                    @test_throws MethodError convert(T, arr)
                end
            end

            @testset "FixedSizeMatrix" begin
                m = FSM{Float64}(undef, 3, 3)
                m_ref = reshape(1:9, size(m))
                @test length(m) == 9
                test_inferred_noalloc(((a, b) -> a .= b), FSM{Float64}, (m, m_ref))
                @test m == m_ref
                test_inferred(similar, FSM{Float64}, (m,))
                test_inferred(similar, FSM{Int}, (m, Int))
                test_inferred(copy, FSM{Float64}, (m,))
                test_inferred(zero, FSM{Float64}, (m,))
                test_inferred(similar, FSM{Int}, (FSM{Int}, (2, 3)))
                test_inferred(similar, FSM{Int}, (FSA{Int}, (2, 3)))
                example_abstract_matrices = (
                    reshape(1:9, (3, 3)),
                    OffsetArray(reshape(1:9, (3, 3)), 1:3, 1:3),
                )
                for T in (FSA{Int}, FSM{Int}), arr in example_abstract_matrices
                    # Don't allow conversion:
                    # <https://github.com/JuliaArrays/FixedSizeArrays.jl/issues/120>.
                    @test_throws MethodError convert(T, arr)
                end
            end

            @testset "`map`" begin
                for dim_count ∈ 0:3
                    size = ntuple(Returns(3), dim_count)
                    a = FSA{Int, dim_count}(undef, size)
                    for v ∈ (3, 3.1, nothing)
                        args = (Returns(v), a)
                        @test all(==(v), map(args...))
                        test_inferred(map, FSA{typeof(v), dim_count}, args)
                    end
                end
            end

            @testset "broadcasting" begin
                TVI = FSV{Int}
                TVF = FSV{Float64}
                v3 = TVI(1:3)
                vf3 = TVF(1:3)
                bp = (x, y) -> x .+ y
                bm = (x, y) -> x .* y
                test_inferred(+,  TVI, (v3, v3))
                test_inferred(bp, TVI, (v3, v3))
                test_inferred(bm, TVI, (v3, v3))
                test_inferred(+,  TVF, (v3, vf3))
                test_inferred(bp, TVF, (v3, vf3))
                test_inferred(bm, TVF, (v3, vf3))
                test_inferred(bp, TVI, (v3,  3))
                test_inferred(bm, TVI, (v3,  3))
                test_inferred(bp, TVF, (v3, .3))
                test_inferred(bm, TVF, (v3, .3))
                @testset "views" begin
                    va = @view TVI(1:3)[1:2]
                    vb = @view TVI(4:6)[2:3]
                    vf2 = TVF(7:8)
                    test_inferred(+,  TVI, (va,  vb))
                    test_inferred(+,  TVF, (va, vf2))
                    test_inferred(+,  TVF, (vf2, vb))
                    test_inferred(bp, TVI, (va, vb))
                    test_inferred(bp, TVF, (va, vf2))
                    test_inferred(bp, TVF, (vf2, vb))
                    test_inferred(-,  TVI, (va, vb))
                    test_inferred(-,  TVF, (va, vf2))
                    test_inferred(-,  TVF, (vf2, vb))
                    test_inferred(bm, TVI, (va, vb))
                    test_inferred(bm, TVF, (va, vf2))
                    test_inferred(bm, TVF, (vf2, vb))
                end
                @testset "matrices" begin
                    TMI = FSM{Int}
                    TMF = FSM{Float64}
                    m33 = TMI(undef, 3, 3)
                    m13 = TMI(undef, 1, 3)
                    m31 = TMI(undef, 3, 1)
                    test_inferred(+,  TMI, (m33, m33))
                    test_inferred(bp, TMI, (m33, m33))
                    test_inferred(bm, TMI, (m33, m33))
                    test_inferred(bp, TMI, (m33, m13))
                    test_inferred(bm, TMI, (m33, m13))
                    test_inferred(bp, TMI, (m33, m31))
                    test_inferred(bm, TMI, (m33, m31))
                    test_inferred(bp, TMI, (m33, v3 ))
                    test_inferred(bm, TMI, (m33, v3 ))
                    test_inferred(bp, TMF, (m33, vf3))
                    test_inferred(bm, TMF, (m33, vf3))
                    test_inferred(bp, TMI, (m33,   3))
                    test_inferred(bm, TMI, (m33,   3))
                    test_inferred(bp, TMF, (m33,  .3))
                    test_inferred(bm, TMF, (m33,  .3))
                    vm22 = @view m33[1:2, 2:3]
                    mf22 = TMF(undef, 2, 2)
                    test_inferred(+,  TMI, (vm22, vm22))
                    test_inferred(+,  TMF, (vm22, mf22))
                    test_inferred(+,  TMF, (mf22, vm22))
                    test_inferred(bp, TMI, (vm22, vm22))
                    test_inferred(bp, TMF, (vm22, mf22))
                    test_inferred(bp, TMF, (mf22, vm22))
                    test_inferred(-,  TMI, (vm22, vm22))
                    test_inferred(-,  TMF, (vm22, mf22))
                    test_inferred(-,  TMF, (mf22, vm22))
                    test_inferred(bm, TMI, (vm22, vm22))
                    test_inferred(bm, TMF, (vm22, mf22))
                    test_inferred(bm, TMF, (mf22, vm22))
                end
                @testset "aliasing" begin
                    v = TVI(1:4)
                    c = @inferred Base.unaliascopy(v)
                    @test c isa TVI
                    @test c == v
                    @test parent(c) !== parent(v)
                    # overlapping source and destination are unaliased with `unaliascopy`
                    @views v[2:4] .= v[1:3]
                    @test v == [1, 1, 2, 3]
                end
                @testset "out-of-place" begin
                    # zero-dimensional broadcasting returns a scalar
                    a0 = FSA{Int}(undef, ())
                    a0[] = 3
                    @test (@inferred a0 .+ 1) === 4
                    # empty arrays
                    e = TVF(undef, 0)
                    @test (@inferred e .+ 1) isa TVF
                    @test isempty(e .+ 1)
                    # a non-concrete element type is widened as needed
                    w = FSV{Any}(undef, 2)
                    w[1] = 1
                    w[2] = 2.0
                    r = w .+ 1
                    @test r isa FSV{Real}
                    @test r == [2, 3.0]
                    # different shapes and views of the same array
                    m = FSM{Int}([1 2; 3 4; 5 6])
                    @test (@inferred m .+ v3) isa FSM{Int}
                    @test m .+ v3 == collect(m) .+ collect(v3)
                    @test v3 .* v3' == collect(v3) .* collect(v3)'
                    @test (@views v3[1:2] .+ v3[2:3]) == [3, 5]
                end
            end

            @testset "`isassigned`" begin
                for dim_count ∈ 0:3
                    for elem_type ∈ (Int, FSM{Nothing})
                        size = ntuple(Returns(3), dim_count)
                        a = FSA{elem_type, dim_count}(undef, size)
                        for i_style ∈ (IndexLinear(), IndexCartesian())
                            for i ∈ eachindex(i_style, a)
                                @test isassigned(a, i) == isbitstype(elem_type)
                                test_inferred_noalloc(isassigned, Bool, (a, i))
                            end
                        end
                    end
                end
                @testset "some assigned" begin
                    a = FSM{FSM{Nothing}}(undef, 3, 3)
                    assigned_inds = ((1, 2), (2, 3), (3, 3))
                    for ij ∈ assigned_inds
                        a[ij...] = FSM{Nothing}(undef, 1, 1)
                    end
                    for ij ∈ Iterators.product(1:3, 1:3)
                        @test isassigned(a, ij...) == (ij ∈ assigned_inds)
                        test_inferred_noalloc(isassigned, Bool, (a, ij...))
                    end
                end
            end

            @testset "`copy`" begin
                m = FixedSizeArray(rand(2, 2))::FixedSizeArray
                @test m == copy(m)
                @test m !== copy(m)
                @test reshape(m, :) == m[:]
                @test m !== m[:]
            end

            @testset "`copyto!`" begin
                for (D, S) ∈ (
                    (Vector, FSV),
                    (FSV, Vector),
                    (FSV, FSV),
                    (
                        if @isdefined Memory
                            ((FSV, Memory), (Memory, FSV))
                        else
                            ()
                        end
                    )...,
                )
                    for E ∈ (Float64, Int)
                        s = S{E}(undef, 5)
                        s .= 1:5
                        d = D{Float64}(undef, 5)
                        test_inferred_noalloc(copyto!, D{Float64}, (d, s))
                        @test copyto!(d, s) == 1:5
                        args = (d, 1, s, 1, length(s))
                        test_inferred_noalloc(copyto!, D{Float64}, args)
                        @test copyto!(args...) == 1:5
                        @test_throws ArgumentError copyto!(d, 1, s, 1, -1)
                    end
                end
            end

            @testset "`copy!`" begin
                for S ∈ (Vector, FSV)
                    s = S{Float64}(1:5)
                    d = FSV{Float64}(undef, 5)
                    test_inferred_noalloc(copy!, FSV{Float64}, (d, s))
                    @test copy!(d, s) == 1:5
                    # issue #194: no `resize!` fallback
                    @test_throws ArgumentError copy!(FSV{Float64}(undef, 3), s)
                    @test_throws "arrays must have the same axes for copy! (consider using `copyto!`), got (Base.OneTo(3),) and (Base.OneTo(5),)" copy!(FSV{Float64}(undef, 3), s)
                end
            end

            @testset "LinearAlgebra" begin
                @testset "$T" for T in (Float16, Float32, Float64)
                    v = randn(T, 8)
                    v_fixed = FSV{T}(v)
                    test_inferred(+, FSV{T}, (v_fixed, v_fixed))
                    v_sum = v_fixed + v_fixed
                    @test v_sum ≈ v + v
                    test_inferred(((x, y) -> x * y'), FSM{T}, (v_fixed, v_fixed))
                    v_mul = v_fixed * v_fixed'
                    @test v_mul ≈ v * v'

                    M = randn(T, 8, 8)
                    M_fixed = FSM{T}(M)
                    test_inferred(+, FSM{T}, (M_fixed, M_fixed))
                    M_sum = M_fixed + M_fixed
                    @test M_sum ≈ M + M
                    test_inferred(*, FSM{T}, (M_fixed, M_fixed))
                    M_mul = M_fixed * M_fixed
                    @test M_mul ≈ M * M
                end
            end

            @testset "`pointer`" begin
                for elem_type ∈ (Int8, Int16, Int32, Int64)
                    for dim_count ∈ 0:4
                        siz = ntuple(Returns(2), dim_count)
                        arr = FSA{elem_type, dim_count}(undef, siz)
                        @test pointer(arr) === pointer(arr, 1) === pointer(arr, 2) - sizeof(elem_type)
                        test_inferred_noalloc(pointer, Ptr{elem_type}, (arr,))
                        test_inferred_noalloc(pointer, Ptr{elem_type}, (arr, 1))
                    end
                end
            end

            @testset "`reshape`" begin
                length_to_shapes = Dict(
                    (0 => ((0,), (0, 0), (0, 1), (1, 0), (1, 0, 0), (0, 0, 1))),
                    (1 => ((), (1,), (1, 1), (1, 1, 1))),
                    (2 => ((2,), (1, 2), (2, 1), (1, 2, 1))),
                    (3 => ((3,), (1, 3), (3, 1), (1, 3, 1))),
                    (4 => ((4,), (1, 4), (4, 1), (2, 2), (1, 2, 2), (2, 1, 2))),
                    (6 => ((6,), (1, 6), (6, 1), (2, 3), (3, 2), (1, 3, 2), (2, 1, 3))),
                )
                for elem_type ∈ (Int, Number, Union{Nothing,Int})
                    for len ∈ keys(length_to_shapes)
                        shapes = length_to_shapes[len]
                        for shape1 ∈ shapes
                            a = FSA{elem_type,length(shape1)}(undef, shape1)
                            @test_throws DimensionMismatch reshape(a, length(a)+1)
                            @test_throws DimensionMismatch reshape(a, length(a)+1, 1)
                            @test_throws DimensionMismatch reshape(a, 1, length(a)+1)
                            for shape2 ∈ shapes
                                @test prod(shape1) === prod(shape2) === len  # meta
                                T = FSA{elem_type,length(shape2)}
                                test_inferred_noalloc(reshape, T, (a, shape2))
                                test_inferred(reshape, T, (a, shape2...))
                                test_inferred_noalloc(((a, s) -> reshape(a, s...)), T, (a, shape2))
                                b = reshape(a, shape2)
                                @test size(b) === shape2
                                @test parent(a) === parent(b)
                                @test a === reshape(b, shape1)
                            end
                        end
                    end
                end
            end

            @testset "`collect_as`" begin
                @testset "difficult requested return type" begin
                    T = FixedSizeVectorDefault{T} where {T <: Float32}
                    iter = 3:7
                    # either return a value of the requested type or throw
                    returns = try
                        collect_as(T, iter)
                        true
                    catch e
                        (e isa ArgumentError) || rethrow()
                        false
                    end
                    returns && @test collect_as(T, iter) isa T
                end
                @testset "empty iterator with inexact `eltype`" begin
                    iterator = Iterators.map((x -> x + 0.3), [])
                    @test collect_as(FSV, iterator; empty_iterator_handler = Returns(Union{})) isa FSV{Union{}}
                    @test collect_as(FSV{Union{}}, iterator) isa FSV{Union{}}
                    @test collect_as(FSV{Float32}, iterator) isa FSV{Float32}
                    @test collect_as(FSV{Any}, iterator) isa FSV{Any}
                    @test_throws ArgumentError collect_as(FSV, iterator)
                end
                @testset "explicit handling of empty iterators" begin
                    @test [] == collect_as(FSV, []; empty_iterator_handler = Returns(Union{}))::FixedSizeArray{Union{}}
                    @test [] == collect_as(FSV, []; empty_iterator_handler = Returns(Bool))::FixedSizeArray{Bool}
                end
                @testset "`Union{}`" begin
                    @test_throws Exception collect_as(Union{}, ())
                    @test_throws Exception collect_as(FixedSizeVector{Union{}, Union{}}, ())
                    @test_throws Exception collect_as(FixedSizeVector{<:Any, Union{}}, ())
                end
                for T ∈ (FSA{Int}, FSV{Int})
                    for iterator ∈ (Iterators.repeated(7), Iterators.cycle(7))
                        @test_throws ArgumentError collect_as(T, iterator)
                    end
                end
                @test_throws ErrorException FSA{Int, -1}  # the dimensions type `NTuple{-1, Int}` can not be constructed
                @test_throws TypeError collect_as(FSA{Int, 3.1}, 7:8)
                for T ∈ (FSA{3}, FSV{3})
                    iterator = (7:8, (7, 8))
                    @test_throws MethodError collect_as(T, iterator)
                end
                @testset "buggy iterator with mismatched `size` and `length" begin
                    for iterator ∈ (Iter((), 0, 7), Iter((3, 2), 5, 7))
                        E = eltype(iterator)
                        dim_count = length(size(iterator))
                        for T ∈ (FSA{E}, FSA{E,dim_count})
                            @test_throws DimensionMismatch collect_as(T, iterator)
                        end
                    end
                end
                iterators = (
                    (), (7,), (7, 8), 7, (7 => 8), Ref(7), fill(7),
                    (i for i ∈ 1:3), ((i + 100*j) for i ∈ 1:3, j ∈ 1:2), Iterators.repeated(7, 2),
                    (i for i ∈ 7:9 if i==8), 7:8, 8:7, map(BigInt, 7:8), Int[], [7], [7 8],
                    Iter((), 1, 7), Iter((3,), 3, 7), Iter((3, 2), 6, 7),
                    (i for i ∈ (false, 0x1, 2) if Bool(2 - 1)),
                    (i + false for i ∈ (false, 0x1, 2) if Bool(2 - 2)),
                    Iterators.filter(<(1), Int[]),
                )
                abstract_array_params(::AbstractArray{T,N}) where {T,N} = (T, N)
                @testset "iterator: $iterator" for iterator ∈ iterators
                    a = collect(iterator)
                    (E, dim_count) = abstract_array_params(a)
                    af = collect(Float64, iterator)
                    @test abstract_array_params(af) == (Float64, dim_count)  # meta
                    iszero(dim_count) ||
                    @test_throws DimensionMismatch collect_as(FSA{E,dim_count+1}, iterator)
                    for T ∈ (FSA{E}, FSA{E,dim_count})
                        test_inferred(collect_as, FSA{E,dim_count}, (T, iterator))
                        fsa = collect_as(T, iterator)
                        @test a == fsa
                        @test first(abstract_array_params(fsa)) <: E
                    end
                    # `collect_as` throws by default when given an empty iterator with imprecise `eltype`
                    (isconcretetype(eltype(iterator)) || (eltype(iterator) <: Union{}) || !isempty(iterator)) &&
                    for T ∈ (FSA, FSA{<:Any,dim_count})
                        @test collect_as(T, iterator) isa FSA{E,dim_count}
                        fsa = collect_as(T, iterator)
                        @test a == fsa
                        @test first(abstract_array_params(fsa)) <: E
                    end
                    for T ∈ (FSA{Float64}, FSA{Float64,dim_count})
                        test_inferred(collect_as, FSA{Float64,dim_count}, (T, iterator))
                        fsa = collect_as(T, iterator)
                        @test af == fsa
                        @test first(abstract_array_params(fsa)) <: Float64
                    end
                    # `collect_as` throws by default when given an empty iterator with imprecise `eltype`
                    (isconcretetype(eltype(iterator)) || (eltype(iterator) <: Union{}) || !isempty(iterator)) &&
                    for T ∈ (FixedSizeArray{<:Any,dim_count},)
                        @test collect_as(T, iterator) isa FixedSizeArray{E,dim_count}
                        fsa = collect_as(T, iterator)
                        @test a == fsa
                        @test first(abstract_array_params(fsa)) <: E
                    end
                end
            end
        end
    end

    include("static_dims.jl")

    @testset "Random" begin
        for T in (Float16, Float32, Float64),
            N in (4, 16, 64, 128, 256),
            seed in (17, 42, 123, 5745192337902587512),
            rng in (Random.default_rng(), Random.MersenneTwister())

            v = FixedSizeArrays.default_underlying_storage_type{T}(undef, N)
            fv = FixedSizeVectorDefault{T}(undef, N)
            rng1 = Random.seed!(rng, seed)
            Random.rand!(rng1, v)
            rng2 = Random.seed!(rng, seed)
            # Make sure `rand!` returns `fv` itself.
            @test Random.rand!(rng2, fv) === fv
            # Make sure rand!(::FixedSizeArrays) generates same stream of
            # numbers as rand!(::Memory)
            @test v == fv
        end
    end

    @testset "Adapt" begin
        # Preserver wrapper types:
        arr = transpose([10 20; 30 14])
        farr = adapt(FixedSizeArray, arr)
        @test farr isa Transpose
        @test parent(farr) isa FixedSizeArray
        # A `FixedSizeArray` with a CPU memory backend is converted as a whole
        a = FixedSizeArray([1.0 2.0; 3.0 4.0])
        @test adapt(Array, a) isa Matrix{Float64}
        @test adapt(Array, a) == a
        # An array with static dimensions is rewrapped, independently of the memory backend
        s = wrap(FixedSizeArrayDims, [1.0, 2.0, 3.0, 4.0], (static(2), 2))
        @test adapt(Array, s) isa FixedSizeArrayDims{Float64, 2, Tuple{StaticInt{2}, Int}, Vector{Float64}}
        @test adapt(Array, s) == s
        # but converted to a `FixedSizeArray` as a whole
        f = adapt(FixedSizeArray, s)
        @test f isa FixedSizeArray{Float64, 2}
        @test f == s
        @test parent(f) !== parent(s)
        @test adapt(FixedSizeArray, transpose(s)) isa Transpose{Float64, <:FixedSizeArray}
    end


    @testset "MPI" begin
        MPI.Initialized() || MPI.Init()
        comm = MPI.COMM_WORLD
        rank = MPI.Comm_rank(comm)
        nranks = MPI.Comm_size(comm)
        root = 0
        @testset "$Mem, $T, $sz" for Mem in storage_types,
                                     T in (Float64, Int, ComplexF64, MPITestStruct),
                                     sz in ((6,), (2, 3), (2, 3, 4))
            N = length(sz)
            FSA = FixedSizeArray{T,N,Mem{T}}
            a = FSA(undef, sz)
            for i in eachindex(a)
                a[i] = mpi_test_value(T, 10 * rank + i)
            end
            a0 = copy(a)

            @testset "Buffer and pointer conversion" begin
                buf = MPI.Buffer(a)
                @test buf isa MPI.Buffer{FSA}
                @test buf.data === a
                @test buf.count == length(a)
                @test buf.datatype == MPI.Datatype(T)
                @test MPI.Buffer_send(a).data === a
                c = Base.cconvert(MPI.MPIPtr, a)
                GC.@preserve c begin
                    @test reinterpret(Ptr{T}, Base.unsafe_convert(MPI.MPIPtr, c)) == pointer(a)
                end
            end

            @testset "point-to-point to self" begin
                b = FSA(undef, sz)
                rreq = MPI.Irecv!(b, comm; source = rank, tag = 1)
                sreq = MPI.Isend(a, comm; dest = rank, tag = 1)
                MPI.Waitall([rreq, sreq])
                @test b == a

                c = FSA(undef, sz)
                sreq = MPI.Isend(a, comm; dest = rank, tag = 2)
                MPI.Recv!(c, comm; source = rank, tag = 2)
                MPI.Wait(sreq)
                @test c == a
            end

            @testset "collectives" begin
                MPI.Bcast!(a, comm; root)
                @test a == (rank == root ? a0 : mpi_test_value.(T, 1:length(a)))
                copyto!(a, a0)

                # gather the whole array, and a copy of a slice of it as in
                # `vec(a[1, .., :])`
                slice = N == 1 ? a[1:2] : vec(a[1, ntuple(Returns(:), N - 1)...])
                for sendbuf in (a, slice)
                    @test sendbuf isa FixedSizeArray
                    n = length(sendbuf)
                    recv = rank == root ?
                           MPI.VBuffer(FixedSizeVector{T,Mem{T}}(undef, nranks * n),
                                       fill(n, nranks)) : nothing
                    MPI.Gatherv!(sendbuf, recv, comm; root)
                    if rank == root
                        @test recv.data isa FixedSizeVector
                        @test recv.data[(rank * n + 1):((rank + 1) * n)] == vec(sendbuf)
                    end
                end

                gathered = FixedSizeVector{T,Mem{T}}(undef, nranks * length(a))
                MPI.Allgather!(a, MPI.UBuffer(gathered, length(a)), comm)
                @test gathered[(rank * length(a) + 1):((rank + 1) * length(a))] == vec(a)

                src = rank == root ?
                      MPI.VBuffer(FixedSizeVector{T,Mem{T}}(mpi_test_value.(T, 1:(nranks * length(a)))),
                                  fill(length(a), nranks)) : nothing
                dst = FSA(undef, sz)
                MPI.Scatterv!(src, dst, comm; root)
                @test vec(dst) == mpi_test_value.(T, (rank * length(a) + 1):((rank + 1) * length(a)))
            end

            # user-defined reduction operators are not supported on all platforms
            T == MPITestStruct || @testset "reductions" begin
                b = FSA(undef, sz)
                MPI.Allreduce!(a, b, +, comm)
                expected = sum(r -> mpi_test_value.(T, 10 * r .+ (1:length(a))), 0:(nranks - 1))
                @test vec(b) == expected
                c = MPI.Allreduce(a, +, comm)
                @test c isa FixedSizeArray{T,N}
                @test vec(c) == expected
                MPI.Allreduce!(a, +, comm)  # in place
                @test vec(a) == expected
                copyto!(a, a0)

                d = FSA(undef, sz)
                MPI.Reduce!(a, rank == root ? d : nothing, +, comm; root)
                rank == root && @test vec(d) == expected
                MPI.Reduce!(a, +, comm; root)  # in place on the root
                rank == root && @test vec(a) == expected
                copyto!(a, a0)
            end

            if N > 1
                @testset "strided views" begin
                    # `FastSubArray` (vector datatype), as in `view(a, 1, .., :)`
                    v_fast = view(a, 1, ntuple(Returns(:), N - 1)...)
                    # non-fast `SubArray` with unit ranges (subarray datatype), for which
                    # `MPI.Buffer` stores the parent `FixedSizeArray`
                    v_slow = view(a, 1:1, 2:size(a, 2), ntuple(Returns(:), N - 2)...)
                    @test MPI.Buffer(v_slow).data === a
                    for v in (v_fast, v_slow)
                        b = Array{T}(undef, size(v))
                        rreq = MPI.Irecv!(b, comm; source = rank, tag = 3)
                        sreq = MPI.Isend(v, comm; dest = rank, tag = 3)
                        MPI.Waitall([rreq, sreq])
                        @test b == v

                        # receive into a strided view of a `FixedSizeArray`
                        c = FSA(undef, sz)
                        fill!(c, mpi_test_value(T, 0))
                        w = view(c, parentindices(v)...)
                        rreq = MPI.Irecv!(w, comm; source = rank, tag = 4)
                        sreq = MPI.Isend(b, comm; dest = rank, tag = 4)
                        MPI.Waitall([rreq, sreq])
                        @test w == v
                        @test count(==(mpi_test_value(T, 0)), c) == length(c) - length(v)
                    end

                    # scatter into a strided view
                    n = length(v_fast)
                    src = rank == root ?
                          MPI.VBuffer(FixedSizeVector{T,Mem{T}}(mpi_test_value.(T, 1:(nranks * n))),
                                      fill(n, nranks)) : nothing
                    c = copy(a)
                    MPI.Scatterv!(src, view(c, 1, ntuple(Returns(:), N - 1)...), comm; root)
                    @test vec(view(c, 1, ntuple(Returns(:), N - 1)...)) ==
                          mpi_test_value.(T, (rank * n + 1):((rank + 1) * n))
                    @test view(c, 2:size(c, 1), ntuple(Returns(:), N - 1)...) ==
                          view(a, 2:size(a, 1), ntuple(Returns(:), N - 1)...)
                end
            end
        end

        @testset "static dimensions" begin
            for Mem ∈ storage_types
                a = FixedSizeArrayDims{Float64, 3, Tuple{StaticInt{2}, StaticInt{3}, Int}, Mem{Float64}}(undef, (static(2), static(3), 4))
                vec(a) .= 10 * rank .+ (1:24)
                buf = MPI.Buffer(a)
                @test buf isa MPI.Buffer{typeof(a)}
                @test buf.data === a
                @test buf.count == 24
                @test buf.datatype == MPI.Datatype(Float64)
                c = Base.cconvert(MPI.MPIPtr, a)
                GC.@preserve c begin
                    @test reinterpret(Ptr{Float64}, Base.unsafe_convert(MPI.MPIPtr, c)) == pointer(a)
                end
                b = similar(a)
                rreq = MPI.Irecv!(b, comm; source = rank, tag = 6)
                sreq = MPI.Isend(a, comm; dest = rank, tag = 6)
                MPI.Waitall([rreq, sreq])
                @test b == a
                MPI.Allreduce!(a, b, +, comm)
                @test vec(b) == sum(r -> 10 * r .+ (1:24), 0:(nranks - 1))
            end
        end

        @testset "conversion is delegated to other storage types" begin
            a = FixedSizeArrays.wrap(FixedSizeArray, MPITestVector(collect(1.0:6.0)), (2, 3))
            @test Base.cconvert(MPI.MPIPtr, a) === Base.cconvert(MPI.MPIPtr, parent(a).v)
            b = FixedSizeArray{Float64}(undef, 2, 3)
            rreq = MPI.Irecv!(b, comm; source = rank, tag = 5)
            sreq = MPI.Isend(a, comm; dest = rank, tag = 5)
            MPI.Waitall([rreq, sreq])
            @test b == a
        end
    end

    @testset "Constant-prop through iterate" begin
        @test Base.return_types() do
            m = FixedSizeArray{Int}(undef, 3)
            m[1] = 1
            m[2] = 2
            m[3] = 3
            tuple(m...)
       end[] == NTuple{3, Int}
    end

    @testset "`size` loads all dimensions at once" begin
        # `size` must stay inferable and foldable, while passing each dimension through an
        # intrinsic, see the comment at its definition (and issue #201).
        for dims in ((3,), (2, 3), (2, 3, 4), (5, 4, 4, 4, 2))
            a = FixedSizeArray{Float64}(undef, dims)
            @test size(a) === dims
            test_inferred(size, typeof(dims), (a,))
            @test Core.Compiler.is_foldable(Base.infer_effects(size, (typeof(a),))) skip=!build_is_production_build
        end
    end

    @testset "getindex with a single colon must return an `AbstractVector`" begin
        for m in 0:3
            @test let a = FixedSizeArray{Float32}(undef, ntuple(Returns(1), m))
                a::AbstractArray{<:Any, m}
                a[:] isa AbstractVector
            end
        end
    end
end
