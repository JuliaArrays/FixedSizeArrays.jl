module FixedSizeArraysAdaptExt

using Adapt
using FixedSizeArrays: FixedSizeArray, FixedSizeArrayDims, new_fixed_size_array, static_dims

Adapt.adapt_storage(fixed::Type{<:FixedSizeArray}, A) = fixed(A)

# A `FixedSizeArray` (all dimensions dynamic) is converted as a whole to another array type,
# as other arrays, e.g. `adapt(CuArray, a)` is a `CuArray`.
Adapt.adapt_structure(to::Type{<:AbstractArray}, a::FixedSizeArray) = Adapt.adapt_storage(to, a)
Adapt.adapt_structure(fixed::Type{<:FixedSizeArray}, a::FixedSizeArray) = Adapt.adapt_storage(fixed, a)

# For all other adaptors, and for arrays with static dimensions, adapt the memory backend and
# wrap it again, as for other array wrappers. Adaptors which are not array types usually have
# storage rules only for the arrays they know (e.g. the adaptors converting the arguments of
# a GPU kernel, which turn a `Vector` or a `CuVector` into an `isbits` device array), so that
# they would leave a `FixedSizeArray` unchanged as a whole. With static dimensions this also
# applies to array types, e.g. `adapt(CuArray, a)` moves the memory backend to the GPU, since
# the static dimensions cannot be represented by other array types.
function Adapt.adapt_structure(to, a::FixedSizeArrayDims)
    new_fixed_size_array(adapt(to, parent(a)), static_dims(a))
end

# Convert to a `FixedSizeArray` as a whole, as for other arrays.
Adapt.adapt_structure(fixed::Type{<:FixedSizeArray}, a::FixedSizeArrayDims) = Adapt.adapt_storage(fixed, a)

end # module FixedSizeArraysAdaptExt
