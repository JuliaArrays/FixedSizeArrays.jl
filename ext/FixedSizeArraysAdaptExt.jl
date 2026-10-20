module FixedSizeArraysAdaptExt

using Adapt
using FixedSizeArrays: FixedSizeArray, FixedSizeArrayDims, new_fixed_size_array, static_dims

Adapt.adapt_storage(fixed::Type{<:FixedSizeArray}, A) = fixed(A)

# A `FixedSizeArray` (all dimensions dynamic) is converted as a whole, as other arrays, e.g.
# `adapt(CuArray, a)` is a `CuArray`.
Adapt.adapt_structure(to, a::FixedSizeArray) = Adapt.adapt_storage(to, a)
Adapt.adapt_structure(fixed::Type{<:FixedSizeArray}, a::FixedSizeArray) = Adapt.adapt_storage(fixed, a)

# With static dimensions, adapt the memory backend and wrap it again, e.g. for converting an
# array backed by a GPU vector for use in a kernel, or for moving it to the GPU with
# `adapt(CuArray, a)`. This is independent of the memory backend, since the static
# dimensions cannot be represented by other array types.
function Adapt.adapt_structure(to, a::FixedSizeArrayDims)
    new_fixed_size_array(adapt(to, parent(a)), static_dims(a))
end

# Convert to a `FixedSizeArray` as a whole, as for other arrays.
Adapt.adapt_structure(fixed::Type{<:FixedSizeArray}, a::FixedSizeArrayDims) = Adapt.adapt_storage(fixed, a)

end # module FixedSizeArraysAdaptExt
