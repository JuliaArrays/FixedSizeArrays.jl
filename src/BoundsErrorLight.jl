"""
    BoundsErrorLight <: Exception

Like `BoundsError`, but doesn't store the array, to prevent clashing with escape analysis.

Stores the type of the array and the shape instead.
"""
struct BoundsErrorLight <: Exception
    type_of_array::DataType
    shape::Tuple{Vararg{Int}}
    index::Int
end

function Base.showerror(io::IO, ex::BoundsErrorLight)
    typ = ex.type_of_array
    shp = ex.shape
    ind = ex.index
    print(io, "BoundsErrorLight: attempt to access array of type ")
    show(io, typ)
    print(io, " and shape ")
    show(io, shp)
    print(io, " at index [")
    show(io, ind)
    print(io, ']')
    nothing
end

# inline everything to help escape analysis

function throw_boundserrorlight(A::AbstractArray, i::Int)
    @inline
    throw(BoundsErrorLight(typeof(A), size(A), i))
end

function check_bounds_light(A::AbstractArray, i::Int)
    @inline
    checkbounds(Bool, A, i) || throw_boundserrorlight(A, i)
    nothing
end
