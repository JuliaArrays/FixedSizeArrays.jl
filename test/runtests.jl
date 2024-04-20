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
end
