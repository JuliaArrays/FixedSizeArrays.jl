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
    end
end
