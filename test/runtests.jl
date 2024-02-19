using Test
using FixedSizeArrays

@testset "FixedSizeArrays" begin
    v = FixedSizeVector{Float64}(undef, 3)
    @test length(v) == 3
    v .= 1:3
    @test v == 1:3

    m = FixedSizeMatrix{Float64}(undef, 3, 3)
    @test length(m) == 9
    m[:] .= 1:9
    @test m[:] == 1:9
end
