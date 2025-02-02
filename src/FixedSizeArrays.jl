module FixedSizeArrays

using CheckedSizeProduct

const default_underlying_storage_type = (@isdefined Memory) ? Memory : Vector

const optional_memory = (@isdefined Memory) ? (Memory,) : ()
const optional_atomic_memory = (@isdefined AtomicMemory) ? (AtomicMemory,) : ()
const optional_generic_memory = (@isdefined GenericMemory) ? (GenericMemory,) : ()

include("FixedSizeArray.jl")

if isdefined(Base, :dataids) && (Base.dataids isa Function)
    # This is an internal, non-public function which is nevertheless needed to
    # get good performance in some cases (e.g. broadcasting):
    # <https://github.com/JuliaArrays/FixedSizeArrays.jl/issues/63>.
    Base.dataids(a::FixedSizeArray) = Base.dataids(a.mem)
end

include("collect_as.jl")

end # module FixedSizeArrays
