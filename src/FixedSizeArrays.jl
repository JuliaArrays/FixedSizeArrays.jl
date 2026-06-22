module FixedSizeArrays

using CheckedSizeProduct: checked_size_product
using Collects: Collect, collect_as
using LightBoundsErrors: checkbounds_lightboundserror
using LightBoundsErrors: LightBoundsError as BoundsErrorLight

if (@isdefined Memory)
    include("refarray.jl")
end

const default_underlying_storage_type = (@isdefined Memory) ? UnsafeRefArray : Vector

const optional_memory = (@isdefined Memory) ? (UnsafeRefArray,) : ()
const optional_atomic_memory = (@isdefined AtomicMemory) ? (AtomicMemory,) : ()
const optional_generic_memory = (@isdefined GenericMemory) ? (UnsafeRefArray,) : ()

include("FixedSizeArray.jl")

if isdefined(Base, :dataids) && (Base.dataids isa Function)
    # This is an internal, non-public function which is nevertheless needed to
    # get good performance in some cases (e.g. broadcasting):
    # <https://github.com/JuliaArrays/FixedSizeArrays.jl/issues/63>.
    Base.dataids(a::FixedSizeArray) = Base.dataids(parent(a))
end

include("collect_as.jl")

end # module FixedSizeArrays
