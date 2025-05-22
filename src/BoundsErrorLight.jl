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

function show_index(io::IO, x)
    if x isa Colon
        print(io, ':')
    elseif x isa Base.OneTo
        print(io, "1:")
        show(io, last(x))
    else
        show(io, x)
    end
    nothing
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
    if (ind isa AbstractRange) || (ind isa AbstractString)
        show(io, ind)
    else
        for (i, x) ∈ enumerate(ind)
            if 1 < i
                print(io, ", ")
            end
            show_index(io, x)
        end
    end
    print(io, ']')
    nothing
end

function throw_boundserrorlight(A::AbstractArray, i::Int)
    @noinline
    throw(BoundsErrorLight(typeof(A), size(A), i))
end

function check_bounds_light(A::AbstractArray, i::Int)
    @inline
    checkbounds(Bool, A, i) || @noinline throw_boundserrorlight(A, i)
    nothing
end
