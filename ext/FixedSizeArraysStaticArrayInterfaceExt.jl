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

function StaticArrayInterface.known_strides(::Type{<:FixedSizeArrayDims{<:Any, <:Any, D}}) where {D}
    known_strides(D)
end

# The strides known from the dimensions of type `D`, with `nothing` for the others.
@generated function known_strides(::Type{D}) where {D <: Tuple}
    ks = known_dims(D)
    isempty(ks) && return ()
    Tuple(accumulate((s, k) -> (s === nothing || k === nothing) ? nothing : s * k, (1, Base.front(ks)...)))
end

# The device of the memory backend, e.g. a GPU memory backend is not a CPU pointer (by
# default, any strided array with `isbits` elements is assumed to be a CPU pointer).
function StaticArrayInterface.device(::Type{A}) where {A <: FixedSizeArrayDims}
    StaticArrayInterface.device(parent_type(A))
end

end # module FixedSizeArraysStaticArrayInterfaceExt
