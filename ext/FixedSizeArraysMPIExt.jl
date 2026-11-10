module FixedSizeArraysMPIExt

using MPI: MPI, MPIPtr
using FixedSizeArrays: FixedSizeArrayDims, underlying_storage

# MPI.jl (v0.20.28+) converts any `AbstractArray` to an `MPIPtr` through
# `Base.cconvert(Ptr{T}, x)` (or the pointer type given by `MPI.mpi_ptr_type`, e.g.
# `CuPtr{T}` for GPU arrays).  Delegate to the memory backend (`underlying_storage`, a
# `DenseVector` or a view of the memory of a `GenericMemoryRef` backend) so that MPI.jl's
# support for the storage type is used: since its first element is the first element of
# the array, the pointer to it is the pointer to the array, and the result of `cconvert`,
# which `ccall` keeps rooted, keeps the memory alive.
Base.cconvert(::Type{MPIPtr}, a::FixedSizeArrayDims) = Base.cconvert(MPIPtr, underlying_storage(a))

# A `FixedSizeArrayDims` is a `DenseArray`, i.e. its elements are contiguous in memory in
# column-major order, so it can be sent as `length(a)` elements of type `eltype(a)`.
function MPI.Buffer(a::FixedSizeArrayDims)
    MPI.Buffer(a, Cint(length(a)), MPI.Datatype(eltype(a)))
end

end # module FixedSizeArraysMPIExt
