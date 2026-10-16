module FixedSizeArraysStaticArrayInterfaceExt

# Static sizes and strides of `FixedSizeArrayDims`s for StaticArrayInterface.jl, as used e.g. by
# StrideArraysCore.jl (`PtrArray`), LoopVectorization.jl and Polyester.jl. All other traits
# (e.g. `contiguous_axis`, `stride_rank`, `dense_dims`) have suitable defaults for
# `DenseArray`s.

using StaticArrayInterface: StaticArrayInterface
using Static: StaticInt
using FixedSizeArrays: FixedSizeArrayDims, known_dims, parent_type, static_dims

StaticArrayInterface.static_size(a::FixedSizeArrayDims) = static_dims(a)

# The stride in a dimension is a `StaticInt` if all the dimensions before it are static.
function StaticArrayInterface.static_strides(a::FixedSizeArrayDims)
    Base.front(accumulate(*, (StaticInt(1), static_dims(a)...)))
end

function StaticArrayInterface.known_size(::Type{<:FixedSizeArrayDims{<:Any, <:Any, D}}) where {D}
    known_dims(D)
end

# The strides known from the static dimensions, with `nothing` for the others. (As for
# `static_strides`, not the fallback of StaticArrayInterface.jl, which gives `(1,)` for a
# 0-dimensional array.)
function StaticArrayInterface.known_strides(::Type{<:FixedSizeArrayDims{<:Any, <:Any, D}}) where {D}
    Base.front(accumulate((s, k) -> (s === nothing || k === nothing) ? nothing : s * k, (1, known_dims(D)...)))
end

# The device of the memory backend, e.g. a GPU memory backend is not a CPU pointer (by
# default, any strided array with `isbits` elements is assumed to be a CPU pointer).
function StaticArrayInterface.device(::Type{A}) where {A <: FixedSizeArrayDims}
    StaticArrayInterface.device(parent_type(A))
end

end # module FixedSizeArraysStaticArrayInterfaceExt
