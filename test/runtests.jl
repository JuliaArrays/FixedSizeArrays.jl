using Test
using FixedSizeArrays
using OffsetArrays: OffsetArray
import Aqua

const checked_dims = FixedSizeArrays.checked_dims

function allocated(f::F, args::Vararg{Any,N}) where {F,N}
    @allocated f(args...)
end

@testset "FixedSizeArrays" begin
    @testset "Aqua.jl" begin
        Aqua.test_all(FixedSizeArrays)
    end

    @testset "Constructors" begin
        @test FixedSizeArray{Float64,0}(undef) isa FixedSizeArray{Float64,0}
        @test FixedSizeArray{Float64,0}(undef, ()) isa FixedSizeArray{Float64,0}
        for offset ∈ (-1, 0, 2, 3)
            ax = offset:(offset + 1)
            oa = OffsetArray([10, 20], ax)
            @test_throws DimensionMismatch FixedSizeArray(oa)
            @test_throws DimensionMismatch FixedSizeVector(oa)
            @test_throws DimensionMismatch FixedSizeArray{Int}(oa)
            @test_throws DimensionMismatch FixedSizeVector{Int}(oa)
        end
        @test_throws ArgumentError FixedSizeArray{Float64,1}(undef, -1)
        @test_throws ArgumentError FixedSizeArray{Float64,1}(undef, (-1,))
        @test_throws ArgumentError FixedSizeArray{Float64,2}(undef, -1, -1)
        @test_throws ArgumentError FixedSizeArray{Float64,2}(undef, typemax(Int), typemax(Int))
        @test_throws ArgumentError FixedSizeArray{Float64,3}(undef, typemax(Int), typemax(Int), 2)
        @test_throws ArgumentError FixedSizeArray{Float64,4}(undef, typemax(Int), typemax(Int), 2, 4)
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
                @test factorial(n) == prod(t) == @inferred checked_dims(t)
                @test iszero(allocated(checked_dims, t))
            end
        end
    end

    @testset "FixedSizeVector" begin
        v = FixedSizeVector{Float64}(undef, 3)
        @test length(v) == 3
        v .= 1:3
        @test v == 1:3
        @test typeof(v) == typeof(similar(v)) == FixedSizeVector{Float64}
        @test similar(v, Int) isa FixedSizeVector{Int}
        @test eltype(similar(v, Int)) == Int
        @test copy(v) isa FixedSizeVector{Float64}
        @test zero(v) isa FixedSizeVector{Float64}
        @test similar(FixedSizeVector{Int}, (2,)) isa FixedSizeVector{Int}
        @test similar(FixedSizeArray{Int}, (2,)) isa FixedSizeVector{Int}
        @test FixedSizeArray{Int}(undef, 2) isa FixedSizeVector{Int}
        example_abstract_vectors = (7:9, [7, 8, 9], OffsetArray([7, 8, 9], 1:3))
        for T ∈ (FixedSizeArray, FixedSizeVector)
            for a ∈ example_abstract_vectors
                @test convert(T, a) isa FixedSizeVector{Int}
                @test convert(T, a) == a
                @test convert(T, convert(T, a)) isa FixedSizeVector{Int}
                @test convert(T, convert(T, a)) == a
            end
        end
        for T ∈ (FixedSizeArray{Int}, FixedSizeVector{Int})
            for S ∈ (Int, Float64)
                for a ∈ map((c -> map(S, c)), example_abstract_vectors)
                    @test convert(T, a) isa FixedSizeVector{Int}
                    @test convert(T, a) == a
                    @test convert(T, convert(T, a)) isa FixedSizeVector{Int}
                    @test convert(T, convert(T, a)) == a
                end
            end
        end
    end

    @testset "FixedSizeMatrix" begin
        m = FixedSizeMatrix{Float64}(undef, 3, 3)
        m_ref = reshape(1:9, size(m))
        @test length(m) == 9
        m .= m_ref
        @test m == m_ref
        @test typeof(m) == typeof(similar(m)) == FixedSizeMatrix{Float64}
        @test similar(m, Int) isa FixedSizeMatrix{Int}
        @test eltype(similar(m, Int)) == Int
        @test copy(m) isa FixedSizeMatrix{Float64}
        @test zero(m) isa FixedSizeMatrix{Float64}
        @test similar(FixedSizeMatrix{Int}, (2, 3)) isa FixedSizeMatrix{Int}
        @test similar(FixedSizeArray{Int}, (2, 3)) isa FixedSizeMatrix{Int}
        @test FixedSizeArray{Int}(undef, 2, 3) isa FixedSizeMatrix{Int}
        example_abstract_matrices = (
            reshape(1:9, (3, 3)),
            OffsetArray(reshape(1:9, (3, 3)), 1:3, 1:3),
        )
        for T ∈ (FixedSizeArray, FixedSizeMatrix)
            for a ∈ example_abstract_matrices
                @test convert(T, a) isa FixedSizeMatrix{Int}
                @test convert(T, a) == a
                @test convert(T, convert(T, a)) isa FixedSizeMatrix{Int}
                @test convert(T, convert(T, a)) == a
            end
        end
        for T ∈ (FixedSizeArray{Int}, FixedSizeMatrix{Int})
            for S ∈ (Int, Float64)
                for a ∈ map((c -> map(S, c)), example_abstract_matrices)
                    @test convert(T, a) isa FixedSizeMatrix{Int}
                    @test convert(T, a) == a
                    @test convert(T, convert(T, a)) isa FixedSizeMatrix{Int}
                    @test convert(T, convert(T, a)) == a
                end
            end
        end
    end

    @testset "`map`" begin
        for dim_count ∈ 0:3
            size = ntuple(Returns(3), dim_count)
            a = FixedSizeArray{Int, dim_count}(undef, size)
            for v ∈ (3, 3.1, nothing)
                mapped = map(Returns(v), a)
                @test all(==(v), mapped)
                @test mapped isa FixedSizeArray{typeof(v), dim_count}
            end
        end
    end

    @testset "broadcasting" begin
        v3 = FixedSizeArray{Int}(undef, 3)
        @test v3 isa FixedSizeVector{Int}
        @test (@inferred (v3  + v3)) isa FixedSizeVector{Int}
        @test (@inferred (v3 .+ v3)) isa FixedSizeVector{Int}
        @test (@inferred (v3 .* v3)) isa FixedSizeVector{Int}
        @test (@inferred (v3 .+  3)) isa FixedSizeVector{Int}
        @test (@inferred (v3 .*  3)) isa FixedSizeVector{Int}
        @test (@inferred (v3 .+ .3)) isa FixedSizeVector{Float64}
        @test (@inferred (v3 .* .3)) isa FixedSizeVector{Float64}
        @testset "matrices" begin
            m33 = FixedSizeArray{Int}(undef, 3, 3)
            m13 = FixedSizeArray{Int}(undef, 1, 3)
            m31 = FixedSizeArray{Int}(undef, 3, 1)
            @test m33 isa FixedSizeMatrix{Int}
            @test m13 isa FixedSizeMatrix{Int}
            @test m31 isa FixedSizeMatrix{Int}
            @test (@inferred (m33 .+ .3 )) isa FixedSizeMatrix{Float64}
            @test (@inferred (m33 .+  3 )) isa FixedSizeMatrix{Int}
            @test (@inferred (m33 .+ v3 )) isa FixedSizeMatrix{Int}
            @test (@inferred (m33 .+ m13)) isa FixedSizeMatrix{Int}
            @test (@inferred (m33 .+ m31)) isa FixedSizeMatrix{Int}
            @test (@inferred (m33 .+ m33)) isa FixedSizeMatrix{Int}
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
                @test copyto!(d, s) isa D{Float64}
                @test copyto!(d, s) == 1:5
                @test copyto!(d, 1, s, 1, length(s)) isa D{Float64}
                @test copyto!(d, 1, s, 1, length(s)) == 1:5
                @test_throws ArgumentError copyto!(d, 1, s, 1, -1)
            end
        end
    end

    @testset "LinearAlgebra" begin
        @testset "$T" for T in (Float16, Float32, Float64)
            v = randn(T, 8)
            v_fixed = FixedSizeVector(v)
            v_sum = @inferred v_fixed + v_fixed
            @test v_sum isa FixedSizeVector{T}
            @test v_sum ≈ v + v
            v_mul = @inferred v_fixed * v_fixed'
            @test v_mul isa FixedSizeMatrix{T}
            @test v_mul ≈ v * v'

            M = randn(T, 8, 8)
            M_fixed = FixedSizeMatrix(M)
            M_sum = @inferred M_fixed + M_fixed
            @test M_sum isa FixedSizeMatrix{T}
            @test M_sum ≈ M + M
            M_mul = @inferred M_fixed * M_fixed
            @test M_mul isa FixedSizeMatrix{T}
            @test M_mul ≈ M * M
        end
    end

    @testset "`pointer`" begin
        for elem_type ∈ (Int8, Int16, Int32, Int64)
            for dim_count ∈ 0:4
                siz = ntuple(Returns(2), dim_count)
                arr = FixedSizeArray{elem_type, dim_count}(undef, siz)
                @test pointer(arr) === pointer(arr, 1) === pointer(arr, 2) - sizeof(elem_type)
                @test (@inferred pointer(arr))    isa Ptr{elem_type}
                @test (@inferred pointer(arr, 1)) isa Ptr{elem_type}
                @test iszero(allocated(pointer, arr))
                @test iszero(allocated(pointer, arr, 1))
            end
        end
    end
end
