using Test
using FixedSizeArrays
using OffsetArrays: OffsetArray
import Aqua

const checked_dims = FixedSizeArrays.checked_dims
const collect_as = FixedSizeArrays.collect_as

# helpers for testing for allocation or suboptimal inference

allocated(f::F,      args::Tuple) where {F} = @allocated f(args...)
allocated(::Type{T}, args::Tuple) where {T} = @allocated T(args...)
test_noalloc(f::F,      args::Tuple) where {F} = @test iszero(allocated(f, args))
test_noalloc(::Type{T}, args::Tuple) where {T} = @test iszero(allocated(T, args))
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

@testset "FixedSizeArrays" begin
    @testset "meta" begin
        @test isempty(detect_ambiguities(Main))  # test for ambiguities in this file
    end

    @testset "Aqua.jl" begin
        Aqua.test_all(FixedSizeArrays)
    end

    @testset "Constructors" begin
        for dim_count ∈ 0:4
            siz = ntuple(Returns(2), dim_count)
            T = FixedSizeArray{Float64}
            R = FixedSizeArray{Float64,dim_count}
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
            @test_throws DimensionMismatch FixedSizeArray(oa)
            @test_throws DimensionMismatch FixedSizeVector(oa)
            @test_throws DimensionMismatch FixedSizeArray{Int}(oa)
            @test_throws DimensionMismatch FixedSizeVector{Int}(oa)
        end
        @testset "taken from Julia's test/core.jl" begin
            # inspired by:
            # https://github.com/JuliaLang/julia/blob/83929ad883c97fa4376618c4559b74f6ada7a3ce/test/core.jl#L7211-L7228
            b = prevpow(2, typemax(Int))
            test_inferred(FixedSizeArray{Int}, FixedSizeArray{Int,3}, (undef, 0, b, b))
            test_inferred(FixedSizeArray{Int}, FixedSizeArray{Int,3}, (undef, b, b, 0))
            test_inferred(FixedSizeArray{Int}, FixedSizeArray{Int,5}, (undef, b, b, 0, b, b))
            @test_throws ArgumentError FixedSizeArray{Int}(undef, b, b)
            @test_throws ArgumentError FixedSizeArray{Int}(undef, 1, b, b)
            @test_throws ArgumentError FixedSizeArray{Int}(undef, 0, -10)
            @test_throws ArgumentError FixedSizeArray{Int}(undef, -10, 0)
            @test_throws ArgumentError FixedSizeArray{Int}(undef, -1, -1)
            @test_throws ArgumentError FixedSizeArray{Int}(undef, 0, -4, -4)
            @test_throws ArgumentError FixedSizeArray{Int}(undef, -4, 1, 0)
            @test_throws ArgumentError FixedSizeArray{Int}(undef, -4, -4, 1)
        end
        @test_throws ArgumentError FixedSizeArray{Float64,1}(undef, -1)
        @test_throws ArgumentError FixedSizeArray{Float64,1}(undef, (-1,))
        @test_throws ArgumentError FixedSizeArray{Float64,2}(undef, -1, -1)
        @test_throws ArgumentError FixedSizeArray{Float64,2}(undef, typemax(Int), 0)
        @test_throws ArgumentError FixedSizeArray{Float64,2}(undef, typemax(Int), typemax(Int))
        @test_throws ArgumentError FixedSizeArray{Float64,3}(undef, typemax(Int)-1, typemax(Int)-1, 2)
        @test_throws ArgumentError FixedSizeArray{Float64,4}(undef, typemax(Int)-1, typemax(Int)-1, 2, 4)
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
                    @test_throws ArgumentError FixedSizeArray{Float64,n}(undef, siz)
                    @test_throws ArgumentError FixedSizeArray{Float64,n}(undef, siz...)
                end
            end
        end
    end

    @testset "safe computation of length from dimensions size" begin
        @test isone(checked_dims(()))
        for n ∈ 0:30
            t = Tuple(1:n)
            if 20 < n
                @test_throws ArgumentError checked_dims(t)
            else
                @test factorial(n) == prod(t) == checked_dims(t)
                test_inferred_noalloc(checked_dims, Int, (t,))
            end
        end
    end

    @testset "FixedSizeVector" begin
        v = FixedSizeVector{Float64}(undef, 3)
        @test length(v) == 3
        test_inferred_noalloc(((a, b) -> a .= b), FixedSizeVector{Float64}, (v, 1:3))
        @test v == 1:3
        test_inferred(similar, FixedSizeVector{Float64}, (v,))
        test_inferred(similar, FixedSizeVector{Int}, (v, Int))
        test_inferred(copy, FixedSizeVector{Float64}, (v,))
        test_inferred(zero, FixedSizeVector{Float64}, (v,))
        test_inferred(similar, FixedSizeVector{Int}, (FixedSizeVector{Int}, (2,)))
        test_inferred(similar, FixedSizeVector{Int}, (FixedSizeArray{Int}, (2,)))
        example_abstract_vectors = (7:9, [7, 8, 9], OffsetArray([7, 8, 9], 1:3))
        test_convert = function(T, a)
            test_inferred(convert, FixedSizeVector{Int}, (T, a))
            c = convert(T, a)
            test_inferred(convert, FixedSizeVector{Int}, (T, c))
            @test c == a
            @test c === convert(T, c)
        end
        for T ∈ (FixedSizeArray, FixedSizeVector)
            for a ∈ example_abstract_vectors
                test_convert(T, a)
            end
        end
        for T ∈ (FixedSizeArray{Int}, FixedSizeVector{Int})
            for S ∈ (Int, Float64)
                for a ∈ map((c -> map(S, c)), example_abstract_vectors)
                    test_convert(T, a)
                end
            end
        end
    end

    @testset "FixedSizeMatrix" begin
        m = FixedSizeMatrix{Float64}(undef, 3, 3)
        m_ref = reshape(1:9, size(m))
        @test length(m) == 9
        test_inferred_noalloc(((a, b) -> a .= b), FixedSizeMatrix{Float64}, (m, m_ref))
        @test m == m_ref
        test_inferred(similar, FixedSizeMatrix{Float64}, (m,))
        test_inferred(similar, FixedSizeMatrix{Int}, (m, Int))
        test_inferred(copy, FixedSizeMatrix{Float64}, (m,))
        test_inferred(zero, FixedSizeMatrix{Float64}, (m,))
        test_inferred(similar, FixedSizeMatrix{Int}, (FixedSizeMatrix{Int}, (2, 3)))
        test_inferred(similar, FixedSizeMatrix{Int}, (FixedSizeArray{Int}, (2, 3)))
        example_abstract_matrices = (
            reshape(1:9, (3, 3)),
            OffsetArray(reshape(1:9, (3, 3)), 1:3, 1:3),
        )
        test_convert = function(T, a)
            test_inferred(convert, FixedSizeMatrix{Int}, (T, a))
            c = convert(T, a)
            test_inferred(convert, FixedSizeMatrix{Int}, (T, c))
            @test c == a
            @test c === convert(T, c)
        end
        for T ∈ (FixedSizeArray, FixedSizeMatrix)
            for a ∈ example_abstract_matrices
                test_convert(T, a)
            end
        end
        for T ∈ (FixedSizeArray{Int}, FixedSizeMatrix{Int})
            for S ∈ (Int, Float64)
                for a ∈ map((c -> map(S, c)), example_abstract_matrices)
                    test_convert(T, a)
                end
            end
        end
    end

    @testset "`map`" begin
        for dim_count ∈ 0:3
            size = ntuple(Returns(3), dim_count)
            a = FixedSizeArray{Int, dim_count}(undef, size)
            for v ∈ (3, 3.1, nothing)
                args = (Returns(v), a)
                @test all(==(v), map(args...))
                test_inferred(map, FixedSizeArray{typeof(v), dim_count}, args)
            end
        end
    end

    @testset "broadcasting" begin
        TVI = FixedSizeVector{Int}
        TVF = FixedSizeVector{Float64}
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
        @testset "matrices" begin
            TMI = FixedSizeMatrix{Int}
            TMF = FixedSizeMatrix{Float64}
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
        end
    end

    @testset "`isassigned`" begin
        for dim_count ∈ 0:3
            for elem_type ∈ (Int, FixedSizeMatrix{Nothing})
                size = ntuple(Returns(3), dim_count)
                a = FixedSizeArray{elem_type, dim_count}(undef, size)
                for i_style ∈ (IndexLinear(), IndexCartesian())
                    for i ∈ eachindex(i_style, a)
                        @test isassigned(a, i) == isbitstype(elem_type)
                        test_inferred_noalloc(isassigned, Bool, (a, i))
                    end
                end
            end
        end
        @testset "some assigned" begin
            a = FixedSizeMatrix{FixedSizeMatrix{Nothing}}(undef, 3, 3)
            assigned_inds = ((1, 2), (2, 3), (3, 3))
            for ij ∈ assigned_inds
                a[ij...] = FixedSizeMatrix{Nothing}(undef, 1, 1)
            end
            for ij ∈ Iterators.product(1:3, 1:3)
                @test isassigned(a, ij...) == (ij ∈ assigned_inds)
                test_inferred_noalloc(isassigned, Bool, (a, ij...))
            end
        end
    end

    @testset "`copyto!`" begin
        for (D, S) ∈ (
            (Memory, FixedSizeVector),
            (FixedSizeVector, Memory),
            (FixedSizeVector, FixedSizeVector),
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

    @testset "LinearAlgebra" begin
        @testset "$T" for T in (Float16, Float32, Float64)
            v = randn(T, 8)
            v_fixed = FixedSizeVector(v)
            test_inferred(+, FixedSizeVector{T}, (v_fixed, v_fixed))
            v_sum = v_fixed + v_fixed
            @test v_sum ≈ v + v
            test_inferred(((x, y) -> x * y'), FixedSizeMatrix{T}, (v_fixed, v_fixed))
            v_mul = v_fixed * v_fixed'
            @test v_mul ≈ v * v'

            M = randn(T, 8, 8)
            M_fixed = FixedSizeMatrix(M)
            test_inferred(+, FixedSizeMatrix{T}, (M_fixed, M_fixed))
            M_sum = M_fixed + M_fixed
            @test M_sum ≈ M + M
            test_inferred(*, FixedSizeMatrix{T}, (M_fixed, M_fixed))
            M_mul = M_fixed * M_fixed
            @test M_mul ≈ M * M
        end
    end

    @testset "`pointer`" begin
        for elem_type ∈ (Int8, Int16, Int32, Int64)
            for dim_count ∈ 0:4
                siz = ntuple(Returns(2), dim_count)
                arr = FixedSizeArray{elem_type, dim_count}(undef, siz)
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
                    a = FixedSizeArray{elem_type,length(shape1)}(undef, shape1)
                    @test_throws DimensionMismatch reshape(a, length(a)+1)
                    @test_throws DimensionMismatch reshape(a, length(a)+1, 1)
                    @test_throws DimensionMismatch reshape(a, 1, length(a)+1)
                    for shape2 ∈ shapes
                        @test prod(shape1) === prod(shape2) === len  # meta
                        T = FixedSizeArray{elem_type,length(shape2)}
                        test_inferred_noalloc(reshape, T, (a, shape2))
                        test_inferred_noalloc(reshape, T, (a, shape2...))
                        b = reshape(a, shape2)
                        @test size(b) === shape2
                        @test a.mem === b.mem
                        @test a === reshape(b, shape1)
                    end
                end
            end
        end
    end

    @testset "`collect_as`" begin
        for T ∈ (FixedSizeArray, FixedSizeVector, FixedSizeArray{Int}, FixedSizeVector{Int})
            for iterator ∈ (Iterators.repeated(7), Iterators.cycle(7))
                @test_throws ArgumentError collect_as(T, iterator)
            end
        end
        for T ∈ (FixedSizeArray{<:Any,-1}, FixedSizeArray{Int,-1}, FixedSizeArray{Int,3.1})
            iterator = (7:8, (7, 8))
            @test_throws ArgumentError collect_as(T, iterator)
        end
        for T ∈ (FixedSizeArray{3}, FixedSizeVector{3})
            iterator = (7:8, (7, 8))
            @test_throws TypeError collect_as(T, iterator)
        end
        iterators = (
            (), (7,), (7, 8), 7, (7 => 8), Ref(7), fill(7),
            (i for i ∈ 1:3), ((i + 100*j) for i ∈ 1:3, j ∈ 1:2), Iterators.repeated(7, 2),
            (i for i ∈ 7:9 if i==8), 7:8, 8:7, Int[], [7], [7 8],
        )
        abstract_array_params(::AbstractArray{T,N}) where {T,N} = (T, N)
        @testset "iterator: $iterator" for iterator ∈ iterators
            a = collect(iterator)
            (E, dim_count) = abstract_array_params(a)
            af = collect(Float64, iterator)
            @test abstract_array_params(af) == (Float64, dim_count)  # meta
            @test_throws MethodError collect_as(FixedSizeArray{E,dim_count+1}, iterator)
            for T ∈ (FixedSizeArray, FixedSizeArray{<:Any,dim_count})
                @test a == collect_as(T, iterator)
                fsa = collect_as(T, iterator)
                @test first(abstract_array_params(fsa)) <: E
            end
            for T ∈ (FixedSizeArray{E}, FixedSizeArray{E,dim_count})
                @test a == collect_as(T, iterator)
                test_inferred(collect_as, FixedSizeArray{E,dim_count}, (T, iterator))
                fsa = collect_as(T, iterator)
                @test first(abstract_array_params(fsa)) <: E
            end
            for T ∈ (FixedSizeArray{Float64}, FixedSizeArray{Float64,dim_count})
                @test af == collect_as(T, iterator)
                test_inferred(collect_as, FixedSizeArray{Float64,dim_count}, (T, iterator))
                fsa = collect_as(T, iterator)
                @test first(abstract_array_params(fsa)) <: Float64
            end
        end
    end
end
