module FixedSizeArraysStaticExt

# Static dimensions of `FixedSizeArrayDims`s: `StaticInt`s from Static.jl.

using Static: StaticInt
using FixedSizeArrays: FixedSizeArrays, dim_value, to_static_dim

FixedSizeArrays.dim_value(::StaticInt{N}) where {N} = N

FixedSizeArrays.known_dim(::Type{StaticInt{N}}) where {N} = N

FixedSizeArrays.to_static_dim(d::StaticInt) = d

function FixedSizeArrays.convert_dim(::Type{StaticInt{N}}, d) where {N}
    if dim_value(to_static_dim(d)) != N  # (`to_static_dim`: `Integer`s which are not `Int`s)
        throw(DimensionMismatch(lazy"the dimension $d does not match the static dimension $N"))
    end
    StaticInt(N)
end

end # module FixedSizeArraysStaticExt
