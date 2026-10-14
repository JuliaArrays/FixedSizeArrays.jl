module FixedSizeArraysAdaptExt

using Adapt
using FixedSizeArrays: FixedSizeArray, StaticSizeArray, new_static_size_array, static_dims

Adapt.adapt_storage(fixed::Type{<:FixedSizeArray}, A) = fixed(A) 

# Adapt the memory backend of a `StaticSizeArray` and wrap it again, e.g. for converting
# a `StaticSizeArray` backed by a GPU vector for use in a kernel, or for moving it to the GPU
# with `adapt(CuArray, a)`. Unlike for `FixedSizeArray`s, this is independent of the memory
# backend, since the static dimensions cannot be represented by other array types.
function Adapt.adapt_structure(to, a::StaticSizeArray)
    new_static_size_array(adapt(to, parent(a)), static_dims(a))
end

# Convert to a `FixedSizeArray` as a whole, as for other arrays.
Adapt.adapt_structure(fixed::Type{<:FixedSizeArray}, a::StaticSizeArray) = Adapt.adapt_storage(fixed, a)

end # module FixedSizeArraysAdaptExt
