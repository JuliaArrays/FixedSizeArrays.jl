export collect_as

function throw_bottom_type()
    throw(ArgumentError("`Union{}` not expected"))
end

function collect_as_storage_type_helper(::Type{Storage}, ::Type) where {S, Storage <: AbstractVector{S}}
    Storage
end
function collect_as_storage_type_helper(::Type{Storage}, ::Type{E}) where {Storage <: AbstractVector, E}
    Storage{E}
end

function make_abstract_vector_from_tuple(::Type{V}, elems::Tuple) where {E, V <: AbstractVector{E}}
    ret = V(undef, length(elems))
    copyto!(ret, elems)
end

function fsv_type_from_underlying_storage_type(::Type{V}) where {E, V <: DenseVector{E}}
    FixedSizeVector{E, V}
end

function make_vector_from_tuple(::Type{V}, elems::Tuple) where {V <: DenseVector}
    stor = collect_as_storage_type_helper(V, eltype(elems))
    ret_type = Vector{eltype(stor)}
    make_abstract_vector_from_tuple(ret_type, elems)
end

function parent_type_with_default(::Type{<:(FixedSizeArray{E, N, T} where {N})}) where {E, T <: DenseVector{E}}
    T
end
function parent_type_with_default(::Type{<:FixedSizeArray{E}}) where {E}
    default_underlying_storage_type{E}
end
function parent_type_with_default(::Type{<:FixedSizeArray})
    default_underlying_storage_type
end
for T ∈ (Vector, optional_memory...)
    FSA = FixedSizeArray{E, N, T{E}} where {E, N}
    @eval begin
        function parent_type_with_default(::Type{$FSA})
            $T
        end
        function parent_type_with_default(::Type{($FSA){E, N} where {E}}) where {N}
            $T
        end
    end
end

function push(v::Vector, e)
    E = typejoin(typeof(e), eltype(v))
    ret = Vector{E}(undef, length(v) + 1)
    ret = copyto!(ret, v)
    ret[end] = e
    ret
end

function push!!(::Type{<:AbstractVector{E}}, v::Vector{E}, e::E) where {E}
    push!(v, e)
end
function push!!(::Type{<:AbstractVector{E}}, v::Vector{E}, e) where {E}
    push!(v, e)
end
function push!!(::Type{<:AbstractVector}, v::Vector{E}, e::E) where {E}
    push!(v, e)
end
function push!!(::Type{<:AbstractVector}, v::Vector, e)
    push(v, e)
end

function empty_fsv(::Type{V}, ::Any) where {E, V <: DenseVector{E}}
    fsv_type_from_underlying_storage_type(V)(undef, 0)
end

function empty_fsv(::Type{V}, iterator) where {V <: DenseVector}
    let E = eltype(iterator)
        if isconcretetype(E)
            fsv_type_from_underlying_storage_type(V{E})(undef, 0)
        else
            fsv_type_from_underlying_storage_type(V{Union{}})(undef, 0)
        end
    end
end

"""
    collect_as_fsv(V::Type{<:DenseVector}, iterator)

Collect the elements of `iterator` into a `FixedSizeVector`. The argument `V`
specifies the underlying storage type parameter of this `FixedSizeVector`. When
possible, the element type of the return value is also taken from `V`, otherwise it
is determined by the types of the elements of the iterator.

When `V` does not provide an element type and `isempty(iterator)`, the element type
of the return value is:
* `eltype(iterator)`, when it's a concrete type
* `Union{}`, otherwise
"""
function collect_as_fsv(::Type{V}, iterator) where {V <: DenseVector}
    es1 = iterate(iterator)  # unroll a bit to avoid unnecessary allocations and help inference
    T2 = Tuple{Any, Any}
    if es1 isa T2
        let (e1, s1) = es1, state = s1, ret = make_vector_from_tuple(V, (e1,))
            while true
                es = iterate(iterator, state)
                if es isa T2
                    let (e, s) = es
                        state = s
                        ret = push!!(V, ret, e)
                    end
                else
                    break
                end
            end
            ret_type_storage = collect_as_storage_type_helper(V, eltype(ret))
            ret_type = fsv_type_from_underlying_storage_type(ret_type_storage)
            ret_type(ret)
        end
    else
        empty_fsv(V, iterator)
    end
end

function checked_dimension_count_of(::Type{<:AbstractArray{<:Any, N}}, input_size_class) where {N}
    n = check_count_value(N)
    m = dimension_count_of(input_size_class)::Int
    if n != m
        throw(DimensionMismatch())
    end
    n
end
function checked_dimension_count_of(::Type, input_size_class)
    check_count_value(dimension_count_of(input_size_class))
end

function collect_as(::Type{<:(FixedSizeArray{E, N, Union{}} where {E, N})}, ::Any)
    throw_bottom_type()
end
function collect_as(::Type{Union{}}, ::Any)
    throw_bottom_type()
end

"""
    collect_as(t::Type{<:FixedSizeArray}, iterator)

Tries to construct a value of type `t` from the iterator `iterator`. The type `t`
must either be concrete, or a `UnionAll` without constraints.
"""
function collect_as(::Type{T}, iterator) where {T<:FixedSizeArray}
    size_class = Base.IteratorSize(iterator)
    if size_class == Base.IsInfinite()
        throw(ArgumentError("iterator is infinite, can't fit infinitely many elements into a `FixedSizeArray`"))
    end
    mem = parent_type_with_default(T)
    output_dimension_count = checked_dimension_count_of(T, size_class)
    fsv = collect_as_fsv(mem, iterator)
    if isone(output_dimension_count)
        fsv
    else
        reshape(fsv, size(iterator))
    end::T
end
