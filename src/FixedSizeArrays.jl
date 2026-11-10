module FixedSizeArrays

using CheckedSizeProduct: checked_size_product
using Collects: Collect, collect_as
using LightBoundsErrors: checkbounds_lightboundserror
using LightBoundsErrors: LightBoundsError
# We use `const` instead of `using...as...` syntax only to work around a bug in
# ExplicitImports: <https://github.com/JuliaTesting/ExplicitImports.jl/issues/185>.
const BoundsErrorLight = LightBoundsError

const default_underlying_storage_type = (@isdefined Memory) ? Memory : Vector

const optional_memory = (@isdefined Memory) ? (Memory,) : ()
const optional_atomic_memory = (@isdefined AtomicMemory) ? (AtomicMemory,) : ()
const optional_generic_memory = (@isdefined GenericMemory) ? (GenericMemory,) : ()
const optional_memory_ref = (@isdefined MemoryRef) ? (MemoryRef,) : ()

# The type of the memory backend (the `Mem` type parameter of `FixedSizeArray{T,N,Mem}`):
# a `DenseVector{T}` (e.g. `Memory{T}` or `Vector{T}`) or, on Julia v1.11 and later, a
# `GenericMemoryRef{kind,T}` (e.g. `MemoryRef{T}`) pointing to the first element.
@static if @isdefined GenericMemoryRef
    const StorageBound{T} = Union{DenseVector{T}, (GenericMemoryRef{K, T, AS} where {K, AS})}
    const StorageBoundAny = Union{DenseVector, GenericMemoryRef}
else
    const StorageBound{T} = DenseVector{T}
    const StorageBoundAny = DenseVector
end

include("FixedSizeArray.jl")

if isdefined(Base, :dataids) && (Base.dataids isa Function)
    # This is an internal, non-public function which is nevertheless needed to
    # get good performance in some cases (e.g. broadcasting):
    # <https://github.com/JuliaArrays/FixedSizeArrays.jl/issues/63>.
    Base.dataids(a::FixedSizeArrayDims) = storage_dataids(parent(a))
    storage_dataids(m) = Base.dataids(m)
end

if isdefined(Base, :unaliascopy) && (Base.unaliascopy isa Function)
    # This is an internal, non-public function.  The generic fallback for
    # `AbstractArray`s is `@noinline`, which leaves an opaque call in every
    # broadcast and penalises e.g. automatic differentiation with Enzyme.
    Base.unaliascopy(a::FixedSizeArray) = copy(a)
end

include("collect_as.jl")

@static if @isdefined GenericMemoryRef
    include("memoryref.jl")
end

@static if VERSION >= v"1.11.0-DEV.469"
    eval(Expr(:public, :wrap, :FixedSizeArrayDims))
end

end # module FixedSizeArrays
