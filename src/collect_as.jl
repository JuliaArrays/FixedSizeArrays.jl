export collect_as

function throw_bottom_type()
    throw(ArgumentError("`Union{}` not expected"))
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

function push_vector(::Type{V}, v, e) where {V <: AbstractVector}
    ret = V(undef, length(v) + 1)
    copyto!(ret, v)
    ret[end] = e
    ret
end

function push_fsv(::Type{V}, v, e) where {E, V <: DenseVector{E}}
    FSV = FixedSizeVector{E, V}
    push_vector(FSV, v, e)
end

function push_fsv(::Type{V}, v, e) where {V <: DenseVector}
    E = typejoin(typeof(e), eltype(v))
    FSV = FixedSizeVector{E, V{E}}
    push_vector(FSV, v, e)
end

function empty_fsv(::Type{V}, ::Any) where {E, V <: DenseVector{E}}
    FixedSizeVector{E, V}(undef, 0)
end

function empty_fsv(::Type{V}, iterator) where {V <: DenseVector}
    let E = eltype(iterator)
        if isconcretetype(E)
            FixedSizeVector{E, V{E}}(undef, 0)
        else
            let E = Union{}
                FixedSizeVector{E, V{E}}(undef, 0)
            end
        end
    end
end

function collect_as_fsv(::Type{V}, iterator) where {V <: DenseVector}
    es1 = iterate(iterator)  # unroll a bit to avoid unnecessary allocations and help inference
    T2 = Tuple{Any, Any}
    if es1 isa T2
        let (e1, s1) = es1, es2 = iterate(iterator, s1)
            if es2 isa T2
                let (e2, s2) = es2, state = s2, ret = push_fsv(V, (e1,), e2)
                    while true
                        es = iterate(iterator, state)
                        if es isa T2
                            let (e, s) = es
                                state = s
                                ret = push_fsv(V, ret, e)
                            end
                        else
                            break
                        end
                    end
                    ret
                end
            else
                push_fsv(V, (), e1)
            end
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
    end
end
