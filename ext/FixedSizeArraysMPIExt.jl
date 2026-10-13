module FixedSizeArraysMPIExt

using MPI: MPI, MPIPtr
using FixedSizeArrays: FixedSizeArray

# MPI.jl (v0.20.28+) converts any `AbstractArray` to an `MPIPtr` through
# `Base.cconvert(Ptr{T}, x)` (or the pointer type given by `MPI.mpi_ptr_type`, e.g.
# `CuPtr{T}` for GPU arrays).  Delegate to the parent so that MPI.jl's support for the
# storage type is used: since a `FixedSizeArray` is always backed by the whole parent
# `DenseVector` starting at its first element, the pointer to the parent is the pointer to
# the array, and the result of `cconvert`, which `ccall` keeps rooted, keeps the parent
# alive.
Base.cconvert(::Type{MPIPtr}, a::FixedSizeArray) = Base.cconvert(MPIPtr, parent(a))

# A `FixedSizeArray` is a `DenseArray`, i.e. its elements are contiguous in memory in
# column-major order, so it can be sent as `length(a)` elements of type `eltype(a)`.
function MPI.Buffer(a::FixedSizeArray)
    MPI.Buffer(a, Cint(length(a)), MPI.Datatype(eltype(a)))
end

end # module FixedSizeArraysMPIExt
