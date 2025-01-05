public collect_as

function collect_as_fsa0(iterator, ::SpecFSA{0,V}) where {V}
    V::UnionAll
    x = only(iterator)
    E = typeof(x)::Type
    ret = FixedSizeArray{E,0,V{E}}(undef)
    ret[] = x
    ret
end

function collect_as_fsa0(iterator, ::SpecFSA{0,V}) where {E,V<:DenseVector{E}}
    E::Type
    x = only(iterator)
    ret = FixedSizeArray{E,0,V}(undef)
    ret[] = x
    ret
end

function fill_fsa_from_iterator!(a, iterator)
    actual_count = 0
    for e ∈ iterator
        actual_count += 1
        a[actual_count] = e
    end
    if actual_count != length(a)
        throw(ArgumentError("`size`-`length` inconsistency"))
    end
end

function collect_as_fsam_with_shape(
    iterator, ::SpecFSA{M,V}, shape::Tuple{Vararg{Int}},
) where {M,V}
    V::UnionAll
    E = eltype(iterator)::Type
    U = V{E}
    ret = FixedSizeArray{E,M,U}(undef, shape)
    fill_fsa_from_iterator!(ret, iterator)
    map(identity, ret)::(FixedSizeArray{T,M,V{T}} where {T})
end

function collect_as_fsam_with_shape(
    iterator, ::SpecFSA{M,V}, shape::Tuple{Vararg{Int}},
) where {M,E,V<:DenseVector{E}}
    E::Type
    ret = FixedSizeArray{E,M,V}(undef, shape)
    fill_fsa_from_iterator!(ret, iterator)
    ret::FixedSizeArray{E,M,V}
end

function collect_as_fsam(iterator, spec::SpecFSA{M}) where {M}
    check_count_value(M)
    shape = if isone(M)
        (length(iterator),)
    else
        size(iterator)
    end::NTuple{M,Any}
    shap = map(Int, shape)::NTuple{M,Int}
    collect_as_fsam_with_shape(iterator, spec, shap)::FixedSizeArray{<:Any,M}
end

function collect_as_fsa1_from_unknown_length(iterator, ::SpecFSA{1,V}) where {V}
    V::UnionAll
    v = collect(iterator)::AbstractVector
    t = FixedSizeVector(v)::FixedSizeVector
    s = map(identity, t)::FixedSizeVector  # fix element type
    et = eltype(s)
    FixedSizeVector{et,V{et}}(s)  # fix underlying storage type
end

function collect_as_fsa1_from_unknown_length(iterator, ::SpecFSA{1,V}) where {E,V<:DenseVector{E}}
    E::Type
    v = collect(E, iterator)::AbstractVector{E}
    T = FixedSizeVector{E,V}
    T(v)::T
end

function collect_as_fsa_impl(iterator, spec::SpecFSA{0}, ::LengthIsKnown)
    collect_as_fsa0(iterator, spec)::FixedSizeArray{<:Any,0}
end

function collect_as_fsa_impl(iterator, spec::SpecFSA, ::LengthIsKnown)
    collect_as_fsam(iterator, spec)::FixedSizeArray
end

function collect_as_fsa_impl(iterator, spec::SpecFSA{1}, ::LengthIsUnknown)
    collect_as_fsa1_from_unknown_length(iterator, spec)::FixedSizeVector
end

function collect_as_fsa_checked(iterator, ::SpecFSA{nothing,V}, ::Val{M}, length_status) where {V,M}
    check_count_value(M)
    collect_as_fsa_impl(iterator, SpecFSA{M,V}(), length_status)::FixedSizeArray{<:Any,M}
end

function collect_as_fsa_checked(iterator, spec::SpecFSA{M}, ::Val{M}, length_status) where {M}
    check_count_value(M)
    collect_as_fsa_impl(iterator, spec, length_status)::FixedSizeArray{<:Any,M}
end

"""
    collect_as(t::Type{<:FixedSizeArray}, iterator)

Tries to construct a value of type `t` from the iterator `iterator`. The type `t`
must either be concrete, or a `UnionAll` without constraints.
"""
function collect_as(::Type{T}, iterator) where {T<:FixedSizeArray}
    spec = fsa_spec_from_type(T)::SpecFSA
    size_class = Base.IteratorSize(iterator)
    if size_class == Base.IsInfinite()
        throw(ArgumentError("iterator is infinite, can't fit infinitely many elements into a `FixedSizeArray`"))
    end
    dim_count_int = dimension_count_of(size_class)
    check_count_value(dim_count_int)
    dim_count = Val(dim_count_int)::Val
    len_stat = length_status(size_class)
    collect_as_fsa_checked(iterator, spec, dim_count, len_stat)::T
end

if isdefined(Base, :dataids) && (Base.dataids isa Function)
    # This is an internal, non-public function which is nevertheless needed to
    # get good performance in some cases (e.g. broadcasting):
    # <https://github.com/JuliaArrays/FixedSizeArrays.jl/issues/63>.
    Base.dataids(a::FixedSizeArray) = Base.dataids(a.mem)
end
