using Test
using FixedSizeArrays

@testset "FixedSizeArrays" begin

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
        for T ∈ (FixedSizeArray, FixedSizeVector)
            a = 1:3
            @test convert(T, a) isa FixedSizeVector{Int}
            @test convert(T, a) == a
            @test convert(T, convert(T, a)) isa FixedSizeVector{Int}
            @test convert(T, convert(T, a)) == a
        end
        for T ∈ (FixedSizeArray{Int}, FixedSizeVector{Int})
            for S ∈ (Int, Float64)
                a = map(S, 1:3)
                @test convert(T, a) isa FixedSizeVector{Int}
                @test convert(T, a) == a
                @test convert(T, convert(T, a)) isa FixedSizeVector{Int}
                @test convert(T, convert(T, a)) == a
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
        for T ∈ (FixedSizeArray, FixedSizeMatrix)
            a = reshape(1:9, (3, 3))
            @test convert(T, a) isa FixedSizeMatrix{Int}
            @test convert(T, a) == a
            @test convert(T, convert(T, a)) isa FixedSizeMatrix{Int}
            @test convert(T, convert(T, a)) == a
        end
        for T ∈ (FixedSizeArray{Int}, FixedSizeMatrix{Int})
            for S ∈ (Int, Float64)
                a = map(S, reshape(1:9, (3, 3)))
                @test convert(T, a) isa FixedSizeMatrix{Int}
                @test convert(T, a) == a
                @test convert(T, convert(T, a)) isa FixedSizeMatrix{Int}
                @test convert(T, convert(T, a)) == a
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
        for dim_count ∈ 1:3
            for elem_type ∈ (Int, FixedSizeMatrix{Nothing})
                size = ntuple(Returns(3), dim_count)
                a = FixedSizeArray{elem_type, dim_count}(undef, size)
                for i_style ∈ (IndexLinear(), IndexCartesian())
                    for i ∈ eachindex(i_style, a)
                        @test isassigned(a, i) == isbitstype(elem_type)
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
    end
end
