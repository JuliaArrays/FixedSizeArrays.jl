module FixedSizeArraysMPIExt

using MPI: MPI, MPIPtr
using FixedSizeArrays: FixedSizeArray, optional_memory

# MPI.jl passes buffers to the C library by `ccall`ing with argument type `MPI.MPIPtr`, so
# that any Julia object `x` with methods for `Base.cconvert(MPIPtr, x)` and
# `Base.unsafe_convert(MPIPtr, cconvert_result)` can be used as a buffer.  MPI.jl defines
# these only for some types (`Array`, `SubArray`, `Ref`, `Ptr`, and GPU arrays in its
# own package extensions), and `MPI.Buffer(x)`, which is called by most high-level
# functions (`Send`, `Recv!`, `Bcast!`, `Gatherv!`, ...), only for `Array`, `Ref` and some
# `SubArray`s.

# Storage types which live in host memory and support `Base.unsafe_convert(Ptr{T}, mem)`.
const HostStorage = Union{Vector, optional_memory...}

# For host memory, return the `FixedSizeArray` itself from `cconvert`: `ccall` keeps the
# result of `cconvert` rooted for the duration of the call, and the `FixedSizeArray` keeps
# its parent alive.  We cannot delegate to `Base.cconvert(MPIPtr, parent(a))` here because
# MPI.jl (up to at least v0.20.27) does not support `Memory`, the default storage type on
# Julia v1.11+.
Base.cconvert(::Type{MPIPtr}, a::FixedSizeArray{T,N,<:HostStorage}) where {T,N} = a

function Base.unsafe_convert(::Type{MPIPtr}, a::FixedSizeArray{T,N,<:HostStorage}) where {T,N}
    reinterpret(MPIPtr, Base.unsafe_convert(Ptr{T}, a))
end

# For any other storage type (e.g. GPU arrays), delegate to the parent so that MPI.jl's
# own support for that array type (if any) is used.  Since a `FixedSizeArray` is always
# backed by the whole parent `DenseVector` starting at its first element, the pointer to
# the parent is the pointer to the array.
Base.cconvert(::Type{MPIPtr}, a::FixedSizeArray) = Base.cconvert(MPIPtr, parent(a))

# A `FixedSizeArray` is a `DenseArray`, i.e. its elements are contiguous in memory in
# column-major order, so it can be sent as `length(a)` elements of type `eltype(a)`.
function MPI.Buffer(a::FixedSizeArray)
    MPI.Buffer(a, Cint(length(a)), MPI.Datatype(eltype(a)))
end

end # module FixedSizeArraysMPIExt
