# Tests of `FixedSizeArrayDims`, the arrays with static dimensions: construction, indexing,
# static sizes and strides, `similar` & co., broadcasting, printing and StaticArrayInterface.jl.
# Included from `runtests.jl`, inside the outer test set, which provides the helpers.

@testset verbose=true "`FixedSizeArrayDims`: static dimensions" begin
    # dimensions `(2, 3, 4)` with all combinations of static and dynamic dimensions
    dims_variants = (
        (static(2), static(3), 4), (static(2), static(3), static(4)), (2, static(3), static(4)),
        (static(2), 3, static(4)), (2, 3, 4), (2, static(3), 4),
    )
    default = FixedSizeArrays.default_underlying_storage_type

    @testset "public" begin
        @test FixedSizeArrayDims === FixedSizeArrays.FixedSizeArrayDims
        @test !Base.isexported(FixedSizeArrays, :FixedSizeArrayDims)
        if VERSION >= v"1.11"
            @test Base.ispublic(FixedSizeArrays, :FixedSizeArrayDims)
        end
    end

    @testset "`FixedSizeArray` is the alias with all dimensions dynamic: $Mem" for Mem ∈ storage_types
        for (T, N) ∈ ((Float64, 0), (Int, 1), (Float32, 2), (Bool, 3))
            R = FixedSizeArray{T, N, Mem{T}}
            @test R === FixedSizeArrayDims{T, N, NTuple{N, Int}, Mem{T}}
            @test isconcretetype(R)
            @test R <: FixedSizeArrayDims{T, N, <:Any, Mem{T}}
            @test !(FixedSizeArrayDims{T, N, <:Any, Mem{T}} <: R)  # (the arrays with static dimensions are not)
        end
        @test FixedSizeArray <: FixedSizeArrayDims
        @test FixedSizeArray{Int} <: FixedSizeArrayDims{Int}
        @test FixedSizeArray{Int, 2} <: FixedSizeArrayDims{Int, 2}
        @test !(FixedSizeArrayDims <: FixedSizeArray)
        @test !(FixedSizeArrayDims{Int, 2, Tuple{StaticInt{2}, Int}, Mem{Int}} <: FixedSizeArray)
        a = FixedSizeArray{Float64}(undef, 2, 3)
        @test a isa FixedSizeArrayDims{Float64, 2, Tuple{Int, Int}, default{Float64}}
        @test FixedSizeArrays.static_dims(a) === (2, 3)
        @test FixedSizeArrays.dim_value.(FixedSizeArrays.static_dims(a)) === (2, 3)
        # the dimensions given to a `FixedSizeArray` constructor are converted to `Int`s
        for (FSA, M) ∈ ((FixedSizeArray{Float64}, default), (FixedSizeArray{Float64, 3}, default), (FixedSizeArray{Float64, 3, Mem{Float64}}, Mem)),
            dims ∈ dims_variants
            b = FSA(undef, dims)
            @test b isa FixedSizeArray{Float64, 3, M{Float64}}
            @test FixedSizeArrays.static_dims(b) === (2, 3, 4)
            test_inferred(FSA, FixedSizeArray{Float64, 3, M{Float64}}, (undef, dims))
        end
        p = Mem{Float64}(undef, 24)
        for dims ∈ (dims_variants..., (static(2), :, 4), (:, 3, static(4)))
            b = wrap(FixedSizeArray, p, dims)
            @test b isa FixedSizeArray{Float64, 3, Mem{Float64}}
            @test FixedSizeArrays.static_dims(b) === (2, 3, 4)
            @test parent(b) === p
            test_inferred(wrap, FixedSizeArray{Float64, 3, Mem{Float64}}, (FixedSizeArray{Float64, 3}, p, dims))
        end
        @test_throws DimensionMismatch FixedSizeArray{Float64, 2}(undef, (static(2), 3, 4))
        @test_throws DimensionMismatch wrap(FixedSizeArray{Float64, 2}, p, (static(2), 3, 4))
        # arrays with static dimensions are not `FixedSizeArray`s
        @test wrap(FixedSizeArrayDims, p, (static(2), 3, 4)) isa FixedSizeArrayDims{Float64, 3, Tuple{StaticInt{2}, Int, Int}, Mem{Float64}}
        @test !(wrap(FixedSizeArrayDims, p, (static(2), 3, 4)) isa FixedSizeArray)
        @test wrap(FixedSizeArrayDims, p, (2, 3, 4)) isa FixedSizeArray{Float64, 3, Mem{Float64}}
    end

    @testset "construction: $Mem" for Mem ∈ storage_types
        p = Mem{Float64}(undef, 24)
        p .= 1:24
        @testset "$dims" for dims ∈ dims_variants
            R = FixedSizeArrayDims{Float64, 3, typeof(dims), Mem{Float64}}
            for SA ∈ (FixedSizeArrayDims, FixedSizeArrayDims{Float64}, FixedSizeArrayDims{Float64, 3}, FixedSizeArrayDims{Float64, 3, D, Mem{Float64}} where {D}, R)
                a = wrap(SA, p, dims)
                @test a isa R
                @test parent(a) === p
                @test size(a) === (2, 3, 4)
                @test FixedSizeArrays.static_dims(a) === dims
                @test a == reshape(1:24, 2, 3, 4)
                test_inferred(wrap, R, (SA, p, dims))
            end
        end
        # `Integer`s are converted to `Int`s
        @test FixedSizeArrays.static_dims(wrap(FixedSizeArrayDims, p, (static(2), Int8(3), UInt(4)))) === (static(2), 3, 4)
        # a `Colon` is computed from the length of the parent
        @test FixedSizeArrays.static_dims(wrap(FixedSizeArrayDims, p, (static(2), :, static(4)))) === (static(2), 3, static(4))
        @test FixedSizeArrays.static_dims(wrap(FixedSizeArrayDims, p, (static(2), static(3), :))) === (static(2), static(3), 4)
        @test_throws DimensionMismatch wrap(FixedSizeArrayDims, p, (static(5), :))
        @test_throws ArgumentError wrap(FixedSizeArrayDims, p, (static(0), :))
        @test_throws ArgumentError wrap(FixedSizeArrayDims, p, (:, :))
        # no copy: writes are visible in both directions
        a = wrap(FixedSizeArrayDims, p, (static(2), static(3), 4))
        a[2, 3, 4] = 100
        @test p[24] == 100
        p[1] = 200
        @test a[1, 1, 1] == 200
        p .= 1:24
        # mismatched sizes
        @test_throws DimensionMismatch wrap(FixedSizeArrayDims, p, (static(2), static(3), 3))
        @test_throws DimensionMismatch wrap(FixedSizeArrayDims, p, (static(2), static(3), 2, 2, 2))
        @test_throws ArgumentError wrap(FixedSizeArrayDims, p, (static(-2), static(3), -4))
        @test_throws ArgumentError wrap(FixedSizeArrayDims, p, (static(2), static(3), -4))
        # type parameters not matching
        @test_throws ArgumentError wrap(FixedSizeArrayDims{Int}, p, (static(2), static(3), 4))
        @test_throws DimensionMismatch wrap(FixedSizeArrayDims{Float64, 2}, p, (static(2), static(3), 4))
        if @isdefined Memory
            other = Mem === Vector ? Memory : Vector
            @test_throws ArgumentError wrap(FixedSizeArrayDims{Float64, 3, D, other{Float64}} where {D}, p, (static(2), static(3), 4))
        end
        # the dimensions are converted to the dimensions type of the array, if it gives it
        S = FixedSizeArrayDims{Float64, 3, Tuple{StaticInt{2}, StaticInt{3}, Int}, Mem{Float64}}
        @test FixedSizeArrays.static_dims(wrap(S, p, (2, static(3), 4))) === (static(2), static(3), 4)
        @test FixedSizeArrays.static_dims(wrap(S, p, (2, 3, static(4)))) === (static(2), static(3), 4)
        @test FixedSizeArrays.static_dims(wrap(S, p, (2, :, 4))) === (static(2), static(3), 4)
        test_inferred(wrap, S, (S, p, (2, 3, 4)))
        @test_throws DimensionMismatch wrap(S, p, (3, static(2), 4))
        @test_throws DimensionMismatch wrap(S, p, (static(3), 2, 4))
        # one dimension
        v = wrap(FixedSizeArrayDims, p, (static(24),))
        @test v isa FixedSizeArrayDims{Float64, 1, Tuple{StaticInt{24}}, Mem{Float64}}
        @test v == 1:24
        # zero dimensions
        z = wrap(FixedSizeArrayDims, Mem{Float64}(undef, 1), ())
        @test z isa FixedSizeArrayDims{Float64, 0, Tuple{}, Mem{Float64}}
        @test z isa FixedSizeArray{Float64, 0, Mem{Float64}}
        @test size(z) === ()
        z[] = 3
        @test z[] == 3
        # zero-length dimensions
        e = wrap(FixedSizeArrayDims, Mem{Float64}(undef, 0), (static(0), static(3), 5))
        @test size(e) == (0, 3, 5)
        @test isempty(e)
        @test size(wrap(FixedSizeArrayDims, Mem{Float64}(undef, 0), (static(2), static(3), 0))) == (2, 3, 0)
        if Mem === Vector
            # re-wrapping a `Vector` after `resize!`
            v = collect(1.0:12.0)
            @test size(wrap(FixedSizeArrayDims, v, (static(2), static(3), :))) == (2, 3, 2)
            resize!(v, 18)
            @test size(wrap(FixedSizeArrayDims, v, (static(2), static(3), :))) == (2, 3, 3)
        end
        # `undef` constructors
        for dims ∈ dims_variants
            D = typeof(dims)
            test_inferred(FixedSizeArrayDims{Int, 3, D, Mem{Int}}, (undef, dims))
            @test size(FixedSizeArrayDims{Int, 3, D, Mem{Int}}(undef, dims)) == (2, 3, 4)
            test_inferred(FixedSizeArrayDims{Int}, FixedSizeArrayDims{Int, 3, D, default{Int}}, (undef, dims))
            test_inferred(FixedSizeArrayDims{Int, 3}, FixedSizeArrayDims{Int, 3, D, default{Int}}, (undef, dims))
            test_inferred(FixedSizeArrayDims{Int, 3, D, Mem{Int}} where {D}, FixedSizeArrayDims{Int, 3, D, Mem{Int}}, (undef, dims))
        end
        @test_throws ArgumentError FixedSizeArrayDims{Int}(undef, (static(2), -1))
        @test_throws DimensionMismatch FixedSizeArrayDims{Int, 2}(undef, (static(2), 3, 4))
        # the static dimensions can be given as `Int`s if the type gives them
        D = Tuple{StaticInt{2}, Int, StaticInt{4}}
        for dims ∈ ((2, 3, 4), (static(2), 3, 4), (2, static(3), static(4)))
            test_inferred(FixedSizeArrayDims{Int, 3, D, Mem{Int}}, (undef, dims))
            # the memory backend is the default one if it is not given
            test_inferred(FixedSizeArrayDims{Int, 3, D}, FixedSizeArrayDims{Int, 3, D, default{Int}}, (undef, dims))
        end
        @test_throws DimensionMismatch FixedSizeArrayDims{Int, 3, D, Mem{Int}}(undef, (3, 3, 4))
        @test_throws DimensionMismatch FixedSizeArrayDims{Int, 2, Tuple{StaticInt{2}, Int}, Mem{Int}}(undef, (2, 3, 4))
        # a partially specified type with a memory backend
        P = FixedSizeArrayDims{Int, 3, D, Mem{Int}} where {D}
        test_inferred(P, FixedSizeArrayDims{Int, 3, Tuple{StaticInt{2}, Int, Int}, Mem{Int}}, (undef, (static(2), 3, 4)))
        # from a `FixedSizeArray`, sharing the memory, and back
        f = wrap(FixedSizeArray, p, (2, 3, 4))
        for dims ∈ dims_variants
            a = wrap(FixedSizeArrayDims, parent(f), dims)
            @test parent(a) === p
            @test a == f
        end
        a = wrap(FixedSizeArrayDims, p, (static(2), static(3), 4))
        b = FixedSizeArray(a)  # copies
        @test b isa FixedSizeArray{Float64, 3}
        @test b == a
        @test parent(b) !== p
        @test wrap(FixedSizeArray, parent(a), size(a)) == a
        @test parent(wrap(FixedSizeArray, parent(a), size(a))) === p
    end

    @testset "indexing: $Mem, $dims" for Mem ∈ storage_types, dims ∈ (
        (static(5), static(4), static(2), 3), (5, static(4), 2, static(3)), (static(5), 4, static(2), 3), (5, 4, 2, 3),
    )
        ref = reshape(collect(1:120), 5, 4, 2, 3)
        p = Mem{Int}(undef, 120)
        p .= vec(ref)
        a = wrap(FixedSizeArrayDims, p, dims)
        @test IndexStyle(a) === IndexLinear()
        @test eachindex(a) === Base.OneTo(120)
        for i ∈ eachindex(ref)
            @test a[i] === ref[i]
        end
        for I ∈ CartesianIndices(ref)
            @test a[I] === ref[I]
            @test a[Tuple(I)...] === ref[I]
            @test a[Tuple(I)..., 1] === ref[I]  # trailing one
        end
        test_inferred_noalloc(getindex, Int, (a, 3, 2, 1, 2))
        test_inferred_noalloc(getindex, Int, (a, 7))
        test_inferred_noalloc(getindex, Int, (a, CartesianIndex(3, 2, 1, 2)))
        b = similar(a)
        for I ∈ CartesianIndices(ref)
            b[Tuple(I)...] = 2 * ref[I]
        end
        @test b == 2 .* ref
        for I ∈ CartesianIndices(ref)
            b[I] = 3 * ref[I]
        end
        @test b == 3 .* ref
        for i ∈ eachindex(ref)
            b[i] = -i
        end
        @test vec(b) == -1:-1:-120
        # out-of-bounds access
        for inds ∈ ((0,), (121,), (6, 1, 1, 1), (1, 5, 1, 1), (1, 1, 3, 1), (1, 1, 1, 4), (0, 1, 1, 1), (1, 1, 1, 1, 2), (1, 1, 1))
            @test_throws BoundsErrorLight a[inds...]
            @test_throws BoundsErrorLight a[inds...] = 1
            @test !checkbounds(Bool, a, inds...)
        end
        @test_throws BoundsErrorLight a[CartesianIndex(1, 1, 1, 4)]
        @test_throws BoundsErrorLight a[CartesianIndex(1, 1, 1, 4)] = 1
        @test_throws ["LightBoundsError: ", "[1, 1, 1, 4]"] a[1, 1, 1, 4]
        # `isassigned`
        @test isassigned(a, 1) && isassigned(a, 5, 4, 2, 3) && !isassigned(a, 121) && !isassigned(a, 1, 1, 1, 4)
        u = wrap(FixedSizeArrayDims, Mem{Any}(undef, 4), (static(2), 2))
        u[2, 1] = 1
        @test [isassigned(u, i, j) for i ∈ 1:2, j ∈ 1:2] == [false false; true false]
        # iteration
        @test collect(a) == ref
        @test [x for x ∈ a] == ref
        @test sum(a) == sum(ref)
    end

    @testset "static sizes and strides" begin
        a = FixedSizeArrayDims{Float64}(undef, (static(5), static(4), static(4), 7))
        T = typeof(a)
        @test size(a) === (5, 4, 4, 7)
        test_inferred_noalloc(size, NTuple{4, Int}, (a,))
        @test Base.return_types(a -> Val(size(a)[1:3]), (T,)) == [Val{(5, 4, 4)}]
        @test Base.return_types(a -> Val(size(a, 2)), (T,)) == [Val{4}]
        @test Base.return_types(a -> Val(size(a, 5)), (T,)) == [Val{1}]
        @test Base.return_types(a -> Val(strides(a)), (T,)) == [Val{(1, 5, 20, 80)}]
        @test Base.return_types(a -> Val(stride(a, 4)), (T,)) == [Val{80}]
        @test Base.return_types(a -> Val(ndims(a)), (T,)) == [Val{4}]
        @test size(a, 4) == 7
        @test size(a, 5) == 1
        @test_throws ArgumentError size(a, 0)
        @test length(a) == 560
        @test strides(a) === (1, 5, 20, 80)
        # static dimensions after a dynamic one
        m = FixedSizeArrayDims{Float64}(undef, (3, static(4), 2, static(5)))
        M = typeof(m)
        @test size(m) === (3, 4, 2, 5)
        @test strides(m) === (1, 3, 12, 24)
        @test Base.return_types(a -> Val(size(a, 2)), (M,)) == [Val{4}]
        @test Base.return_types(a -> Val(size(a, 4)), (M,)) == [Val{5}]
        @test Base.return_types(a -> Val(stride(a, 2)), (M,)) == [Val]
        # the static dimensions take no space
        @test sizeof(T) == sizeof(FixedSizeVectorDefault{Float64})
        @test sizeof(M) == sizeof(FixedSizeArrayDefault{Float64, 2})
        @test Base.elsize(T) == Base.elsize(a) == sizeof(Float64)
        @test pointer(a) === pointer(parent(a))
        @test pointer(a, 2) === pointer(parent(a)) + sizeof(Float64)
        test_inferred_noalloc(pointer, Ptr{Float64}, (a,))
        @test Base.unsafe_convert(Ptr{Float64}, a) === pointer(parent(a))
        @test Base.cconvert(Ptr{Float64}, a) === parent(a)
        @test propertynames(a) === ()
        @test parentindices(a) === (Base.OneTo(560),)
        @test Base.mightalias(a, parent(a))
        @test Base.mightalias(a, wrap(FixedSizeArrayDims, parent(a), (static(5), :)))
        @test !Base.mightalias(a, copy(a))
        # all dimensions dynamic: `size(a, d)` with an `Int` dimension
        f = FixedSizeArray{Float64}(undef, 3, 4)
        @test size(f, 1) == 3 && size(f, 2) == 4 && size(f, 3) == 1
        @test_throws ArgumentError size(f, 0)
        @test Base.return_types(a -> Val(size(a, 3)), (typeof(f),)) == [Val{1}]
        test_inferred_noalloc(size, Int, (f, 1))
    end

    @testset "`similar`, `copy`, `fill!`, `copyto!`, `reshape`: $Mem" for Mem ∈ storage_types
        p = Mem{Float64}(undef, 24)
        p .= 1:24
        dims = (static(2), 3, static(4))
        a = wrap(FixedSizeArrayDims, p, dims)
        R{T} = FixedSizeArrayDims{T, 3, typeof(dims), Mem{T}}
        test_inferred(similar, R{Float64}, (a,))
        test_inferred(similar, R{Int}, (a, Int))
        # with a size: a `FixedSizeArray`, independently of the value of the size
        test_inferred(similar, FixedSizeArray{Int, 3, Mem{Int}}, (a, Int, (2, 3, 4)))
        test_inferred(similar, FixedSizeArray{Float64, 2, Mem{Float64}}, (a, (6, 4)))
        @test size(similar(a, Int, (3, 2, 4))) == (3, 2, 4)
        # range indexing and reductions infer
        test_inferred(a -> a[:, 2:3, :], FixedSizeArray{Float64, 3, Mem{Float64}}, (a,))
        test_inferred(a -> sum(a; dims = 2), FixedSizeArray{Float64, 3, Mem{Float64}}, (a,))
        @test a[:, 2:3, :] == reshape(1:24, 2, 3, 4)[:, 2:3, :]
        @test sum(a; dims = 2) == sum(reshape(1:24, 2, 3, 4); dims = 2)
        # with the type: the static dimensions are kept, from plain `Int`s
        test_inferred(similar, R{Float64}, (R{Float64}, (2, 7, 4)))
        @test size(similar(R{Float64}, (2, 7, 4))) == (2, 7, 4)
        @test_throws DimensionMismatch similar(R{Float64}, (3, 3, 4))
        test_inferred(copy, R{Float64}, (a,))
        @test copy(a) == a
        @test parent(copy(a)) !== p
        test_inferred(zero, R{Float64}, (a,))
        @test iszero(zero(a))
        @test a[:] isa FixedSizeVector{Float64, Mem{Float64}}
        @test a[:] == 1:24
        b = similar(a)
        test_inferred_noalloc(fill!, R{Float64}, (b, 7))
        @test fill!(b, 7) === b
        @test all(==(7), b)
        # `copyto!` with arrays of the same length
        for src ∈ (a, wrap(FixedSizeArray, copy(p), (2, 3, 4)), reshape(collect(1.0:24.0), 2, 3, 4),
                   ((@isdefined Memory) ? (Memory{Float64}(collect(1.0:24.0)),) : ())...)
            fill!(b, 0)
            test_inferred_noalloc(copyto!, R{Float64}, (b, src))
            @test b == a
            fill!(b, 0)
            @test copyto!(b, 1, src, 1, 24) === b
            @test b == a
        end
        for dst ∈ (similar(a), FixedSizeArray{Float64}(undef, 2, 3, 4), zeros(2, 3, 4),
                   ((@isdefined Memory) ? (Memory{Float64}(undef, 24),) : ())...)
            @test copyto!(dst, a) === dst
            @test vec(dst) == 1:24
            @test copyto!(dst, 1, a, 1, 24) === dst
        end
        # `reshape` and `vec` share the memory
        for (r, sz) ∈ ((vec(a), (24,)), (reshape(a, 6, 4), (6, 4)), (reshape(a, (2, 12)), (2, 12)), (reshape(a, :, 2), (12, 2)))
            @test r isa FixedSizeArray{Float64, length(sz), Mem{Float64}}
            @test size(r) == sz
            @test parent(r) === p
        end
        @test_throws DimensionMismatch reshape(a, 5, 5)
    end

    @testset "broadcasting: $Mem" for Mem ∈ storage_types
        ref = reshape(collect(1.0:24.0), 2, 3, 4)
        p = Mem{Float64}(undef, 24)
        p .= vec(ref)
        a = wrap(FixedSizeArrayDims, p, (static(2), static(3), 4))
        F{T} = FixedSizeArray{T, 3, Mem{T}}
        bp = (x, y) -> x .+ y
        test_inferred(bp, F{Float64}, (a, a))
        test_inferred(bp, F{Float64}, (a, 1))
        test_inferred(bp, F{Float64}, (a, ref))
        test_inferred(bp, F{Float64}, (a, wrap(FixedSizeArray, copy(p), (2, 3, 4))))
        test_inferred(bp, F{Float64}, (a, [1.0, 2.0]))
        @test a .+ a == 2 .* ref
        @test a .* [1.0, 2.0] == ref .* [1.0, 2.0]
        b = similar(a)
        test_inferred_noalloc(((x, y) -> x .= y .+ 1), typeof(b), (b, a))
        @test b == ref .+ 1
        b .= 0
        @test iszero(b)
        @test (a .> 12) isa F{Bool}
        @test map(x -> 2x, a) == 2 .* ref
        # views
        v = view(a, :, 2:3, :)
        test_inferred(bp, F{Float64}, (v, v))
        @test v .+ v == 2 .* ref[:, 2:3, :]
        @test Broadcast.BroadcastStyle(typeof(v)) === Broadcast.BroadcastStyle(typeof(a))
        # with different static dimensions
        c = wrap(FixedSizeArrayDims, copy(p), (2, static(3), static(4)))
        test_inferred(bp, F{Float64}, (a, c))
        @test a .+ c == 2 .* ref
    end

    @testset "show" begin
        a = wrap(FixedSizeArrayDims, [1, 2, 3, 4], (static(2), 2))
        @test repr(a) == "[1 3; 2 4]"
        @test sprint(show, MIME"text/plain"(), a) == "2×2 FixedSizeArrayDims{$Int, 2, Tuple{StaticInt{2}, $Int}, Vector{$Int}}:\n 1  3\n 2  4"
        @test sprint(show, typeof(a)) == "FixedSizeArrayDims{$Int, 2, Tuple{StaticInt{2}, $Int}, Vector{$Int}}"
        # A `FixedSizeArray` type is printed through one of the aliases, never as a
        # `FixedSizeArrayDims`, and the printed form round-trips. Which alias is used depends
        # on the Julia version: up to v1.13 Julia prints through an alias only if exactly one
        # matches, and the `make_typealias` method then selects `FixedSizeArray`; later
        # versions select the alias with the fewest parameters, e.g. `FixedSizeMatrixDefault`.
        # (Partially applied types with free trailing parameters, e.g. `FixedSizeArray{Float64}`,
        # are printed with their `where` clause, as for any alias with a bound depending on `T`.)
        for Mem ∈ storage_types
            f = wrap(FixedSizeArray, Mem{Int}(undef, 4), (2, 2))
            for T ∈ (typeof(f), FixedSizeArray{Int, 3, Mem{Int}}, FixedSizeArray{Int, 0, Mem{Int}},
                     Vector{FixedSizeVector{Int, Mem{Int}}}, FixedSizeArray{T, 2, Mem{T}} where {T},
                     FixedSizeArray{T, N, Mem{T}} where {T, N})
                str = sprint(show, T)
                @test !occursin("FixedSizeArrayDims", str)
                @test startswith(str, "FixedSize") || startswith(str, "Vector{FixedSize")
                @test Base.eval(@__MODULE__, Meta.parse(str)) == T
            end
            str = sprint(show, MIME"text/plain"(), f)
            @test !occursin("FixedSizeArrayDims", str)
            @test startswith(str, "2×2 FixedSize")
            @test endswith(str, ":\n" * sprint(Base.print_array, f))
            # without `:compact`, the `struct` is printed, as for `Vector` and `Array`
            @test sprint(show, typeof(f); context = :compact => false) == sprint(show, FixedSizeArrayDims{Int, 2, Tuple{Int, Int}, Mem{Int}}; context = :compact => false)
            @test occursin("FixedSizeArrayDims", sprint(show, typeof(f); context = :compact => false))
        end
        @test sprint(show, FixedSizeArrayDims) == "FixedSizeArrayDims"
        @test sprint(show, FixedSizeArrayDefault) == "FixedSizeArrayDefault"
        @test sprint(show, Union{}) == "Union{}"
    end

    @testset "StaticArrayInterface" begin
        a = FixedSizeArrayDims{Float64}(undef, (static(5), static(4), static(4), 7))
        T = typeof(a)
        @test StaticArrayInterface.static_size(a) === (StaticInt(5), StaticInt(4), StaticInt(4), 7)
        @test StaticArrayInterface.static_strides(a) === (StaticInt(1), StaticInt(5), StaticInt(20), StaticInt(80))
        test_inferred_noalloc(StaticArrayInterface.static_size, Tuple{StaticInt{5}, StaticInt{4}, StaticInt{4}, Int}, (a,))
        test_inferred_noalloc(StaticArrayInterface.static_strides, Tuple{StaticInt{1}, StaticInt{5}, StaticInt{20}, StaticInt{80}}, (a,))
        @test StaticArrayInterface.known_size(T) === (5, 4, 4, nothing)
        @test StaticArrayInterface.known_strides(T) === (1, 5, 20, 80)
        @test StaticArrayInterface.static_size(a, 2) === StaticInt(4)
        @test StaticArrayInterface.static_size(a, 4) === 7
        @test StaticArrayInterface.contiguous_axis(T) === StaticInt(1)
        @test StaticArrayInterface.dense_dims(T) === ntuple(Returns(StaticArrayInterface.True()), 4)
        @test StaticArrayInterface.device(T) === StaticArrayInterface.device(parent(a))
        @test StaticArrayInterface.device(T) === StaticArrayInterface.CPUPointer()
        # static dimensions after a dynamic one
        m = FixedSizeArrayDims{Float64}(undef, (3, static(4), 2, static(5)))
        @test StaticArrayInterface.static_size(m) === (3, static(4), 2, static(5))
        @test StaticArrayInterface.static_strides(m) === (static(1), 3, 12, 24)
        @test StaticArrayInterface.known_size(typeof(m)) === (nothing, 4, nothing, 5)
        @test StaticArrayInterface.known_strides(typeof(m)) === (1, nothing, nothing, nothing)
        v = FixedSizeArrayDims{Float64}(undef, (3,))
        @test v isa FixedSizeVector
        @test StaticArrayInterface.static_size(v) === (3,)
        @test StaticArrayInterface.static_strides(v) === (StaticInt(1),)
        @test StaticArrayInterface.known_strides(typeof(v)) === (1,)
        z = FixedSizeArrayDims{Float64}(undef, ())
        @test StaticArrayInterface.known_size(typeof(z)) === ()
        @test StaticArrayInterface.known_strides(typeof(z)) === ()
        @test StaticArrayInterface.static_strides(z) === ()
        # the device of a `FixedSizeArray` is the one of its memory backend
        @test StaticArrayInterface.device(FixedSizeArray{Float64}(undef, 2)) === StaticArrayInterface.CPUPointer()
        @test StaticArrayInterface.known_size(FixedSizeArray{Float64, 2, Vector{Float64}}) === (nothing, nothing)
    end
end
