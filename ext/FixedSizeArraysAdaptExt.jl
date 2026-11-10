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
    new_fixed_size_array(adapt_backend(to, parent(a), a), static_dims(a))
end

adapt_backend(to, mem, ::FixedSizeArrayDims) = adapt(to, mem)
@static if @isdefined GenericMemoryRef
    # A memory reference is adapted as a `Vector` sharing its elements (without copying). The
    # reference is kept if adapting does not change this `Vector` (e.g. for
    # `adapt(Array, a)`), otherwise the adapted vector (e.g. a GPU vector) is the new backend.
    function adapt_backend(to, ref::MemoryRef, a::FixedSizeArrayDims)
        v = Base.wrap(Array, ref, length(a))
        adapted = adapt(to, v)
        adapted === v ? ref : adapted
    end
end

# Convert to a `FixedSizeArray` as a whole, as for other arrays.
Adapt.adapt_structure(fixed::Type{<:FixedSizeArray}, a::FixedSizeArrayDims) = Adapt.adapt_storage(fixed, a)

end # module FixedSizeArraysAdaptExt
