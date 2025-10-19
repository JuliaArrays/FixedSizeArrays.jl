module FixedSizeArraysAdaptExt

using Adapt
using FixedSizeArrays: FixedSizeArray

Adapt.adapt_storage(fixed::Type{<:FixedSizeArray}, A) = fixed(A) 

end # module FixedSizeArraysAdaptExt
