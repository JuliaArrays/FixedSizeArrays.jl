module FixedSizeArrays

using CheckedSizeProduct: checked_size_product
using Collects: Collect, collect_as
using LightBoundsErrors: checkbounds_lightboundserror
using LightBoundsErrors: LightBoundsError as BoundsErrorLight

const default_underlying_storage_type = (@isdefined Memory) ? Memory : Vector

const optional_memory = (@isdefined Memory) ? (Memory,) : ()
const optional_atomic_memory = (@isdefined AtomicMemory) ? (AtomicMemory,) : ()
const optional_generic_memory = (@isdefined GenericMemory) ? (GenericMemory,) : ()

include("FixedSizeArray.jl")

if isdefined(Base, :dataids) && (Base.dataids isa Function)
    # This is an internal, non-public function which is nevertheless needed to
    # get good performance in some cases (e.g. broadcasting):
    # <https://github.com/JuliaArrays/FixedSizeArrays.jl/issues/63>.
    Base.dataids(a::FixedSizeArrayDims) = Base.dataids(parent(a))
end

if isdefined(Base, :unaliascopy) && (Base.unaliascopy isa Function)
    # This is an internal, non-public function.  The generic fallback for
    # `AbstractArray`s is `@noinline`, which leaves an opaque call in every
    # broadcast and penalises e.g. automatic differentiation with Enzyme.
    Base.unaliascopy(a::FixedSizeArray) = copy(a)
end

include("collect_as.jl")

@static if VERSION >= v"1.11.0-DEV.469"
    eval(Expr(:public, :wrap, :FixedSizeArrayDims))
end

end # module FixedSizeArrays
