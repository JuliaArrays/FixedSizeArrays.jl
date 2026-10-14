module FixedSizeArraysStaticArrayInterfaceExt

# Static sizes and strides of `StaticSizeArray`s for StaticArrayInterface.jl, as used e.g. by
# StrideArraysCore.jl (`PtrArray`), LoopVectorization.jl and Polyester.jl. All other traits
# (e.g. `contiguous_axis`, `stride_rank`, `dense_dims`) have suitable defaults for
# `DenseArray`s.

using StaticArrayInterface: StaticArrayInterface
using Static: known
using FixedSizeArrays: FixedOrStaticSizeArray, StaticSizeArray, parent_type, static_dims, static_strides

StaticArrayInterface.static_size(a::StaticSizeArray) = static_dims(a)

StaticArrayInterface.static_strides(a::StaticSizeArray) = static_strides(static_dims(a))

function StaticArrayInterface.known_size(::Type{<:StaticSizeArray{<:Any, <:Any, D}}) where {D}
    known(D)
end

function StaticArrayInterface.known_strides(::Type{<:StaticSizeArray{<:Any, <:Any, D}}) where {D}
    known_strides(D)
end

# The strides known from the dimensions of type `D`, with `nothing` for the others.
@generated function known_strides(::Type{D}) where {D <: Tuple}
    ks = known(D)
    isempty(ks) && return ()
    Tuple(accumulate((s, k) -> (s === nothing || k === nothing) ? nothing : s * k, (1, Base.front(ks)...)))
end

# The device of the memory backend, e.g. a GPU memory backend is not a CPU pointer (by
# default, any strided array with `isbits` elements is assumed to be a CPU pointer).
function StaticArrayInterface.device(::Type{A}) where {A <: FixedOrStaticSizeArray}
    StaticArrayInterface.device(parent_type(A))
end

end # module FixedSizeArraysStaticArrayInterfaceExt
