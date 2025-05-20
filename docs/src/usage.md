# Usage

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
