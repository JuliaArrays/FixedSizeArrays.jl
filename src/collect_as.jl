function CollectAs.collect_as(::Type{<:(FixedSizeArray{E, N, Union{}} where {E, N})}, ::Any)
    throw(ArgumentError("`Union{}` not expected"))
end

function infer_ndims_impl(::Base.HasShape{N}) where {N}
    N::Int
end

function infer_ndims_impl(::Base.IteratorSize)
    1
end

function infer_ndims(iterator)
    infer_ndims_impl(Base.IteratorSize(iterator))
end

function inferred_shape_impl(output_ndims::Int, iterator)
    if (!isone(output_ndims)) && (output_ndims != infer_ndims(iterator))
        throw(DimensionMismatch("mismatched dimensionalities"))
    end
    if iszero(output_ndims)
        ()
    elseif isone(output_ndims)
        (:,)
    else
        size(iterator)
    end
end

function inferred_shape(::Type{<:(AbstractArray{T, N} where {T})}, iterator) where {N}
    inferred_shape_impl(N::Int, iterator)
end

function inferred_shape(::Type{<:AbstractArray}, iterator)
    inferred_shape_impl(infer_ndims(iterator), iterator)
end

function CollectAs.collect_as(::Type{FSA}, iterator) where {FSA<:FixedSizeArray}
    T = check_constructor_is_allowed(FSA)
    mem = parent_type_with_default(T)
    shape = inferred_shape(T, iterator)
    backing = collect_as(mem, iterator)
    fsv = new_fixed_size_array(backing, (length(backing),))
    reshape(fsv, shape)
end
