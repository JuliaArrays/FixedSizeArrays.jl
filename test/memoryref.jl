# Tests of the `GenericMemoryRef` memory backend (Julia v1.11 and later).

(@isdefined MemoryRef) && @testset "`GenericMemoryRef` backend" begin
    MemRefFSA = FixedSizeArray{T,N,MemoryRef{T}} where {T,N}
    @testset "`wrap` of a reference" begin
        for T in (Int, Float32, String)
            v = T === String ? string.(1:24) : collect(T, 1:24)
            for size in ((24,), (4, 6), (2, 3, 4), (2, 3, 2, 2))
                N = length(size)
                R = FixedSizeArray{T,N,MemoryRef{T}}
                for FSA_ in (FixedSizeArray, FixedSizeArray{T}, FixedSizeArray{T,N}, MemRefFSA, MemRefFSA{T}, R)
                    a = wrap(FSA_, v.ref, size)
                    @test a isa R
                    @test parent(a) === v.ref
                    @test Base.size(a) === size
                    @test a == reshape(v, size)
                    @test wrap(FSA_, v.ref, size...) === a
                    test_inferred(wrap, R, (FSA_, v.ref, size))
                    # `wrap` of the `Vector` with the type asking for the reference
                    if FSA_ <: MemRefFSA
                        @test wrap(FSA_, v, size) === a
                        test_inferred(wrap, R, (FSA_, v, size))
                    end
                end
                @test parentindices(wrap(FixedSizeArray, v.ref, size)) === (Base.OneTo(24),)
            end
            # the `Vector` itself is wrapped unless the reference is asked for
            @test parent(wrap(FixedSizeArray, v, (4, 6))) === v
            @test wrap(MemRefFSA, v) isa FixedSizeVector{T,MemoryRef{T}}
            # writes are visible in both directions
            a = wrap(FixedSizeArray, v.ref, (4, 6))
            a[end] = v[1]
            @test v[end] == v[1]
            v[2] = v[3]
            @test a[2] == v[3]
            @test pointer(a) == pointer(v)
            @test Base.mightalias(a, v)
            @test Base.dataids(a) == Base.dataids(v)
            @test !Base.mightalias(a, copy(a))
            # errors
            @test_throws DimensionMismatch wrap(FixedSizeArray, v.ref, (5, 5))
            @test_throws DimensionMismatch wrap(MemRefFSA, v, 23)  # `Vector`: exact length
            @test_throws ArgumentError wrap(FixedSizeArray, v.ref, (-4, -6))
            @test_throws ArgumentError wrap(FixedSizeArray{Complex{Int}}, v.ref, (4, 6))
            @test_throws DimensionMismatch wrap(FixedSizeArray{T,3}, v.ref, (4, 6))
            @test_throws ArgumentError wrap(FixedSizeArray{T,2,Vector{T}}, v.ref, (4, 6))
            @test_throws ArgumentError wrap(FixedSizeArray{T,2,Memory{T}}, v, (4, 6))
        end
        v = collect(1.0:24.0)
        test_noalloc(p -> wrap(FixedSizeArray, p.ref, (4, 6)), (v,))
        test_noalloc(p -> wrap(FixedSizeArray{Float64,2,MemoryRef{Float64}}, p, (4, 6)), (v,))
    end
    @testset "fewer elements than the memory" begin
        m = Memory{Int}(undef, 10)
        m .= 1:10
        a = wrap(FixedSizeArray, memoryref(m, 3), (2, 3))
        @test a == reshape(3:8, 2, 3)
        @test_throws BoundsErrorLight a[7]
        @test_throws BoundsErrorLight a[0]
        @test_throws BoundsErrorLight a[3, 1]
        @test_throws BoundsErrorLight (a[7] = 1)
        @test !isassigned(a, 7)
        @test wrap(FixedSizeArray, memoryref(m, 3), 8) == 3:10
        @test_throws DimensionMismatch wrap(FixedSizeArray, memoryref(m, 3), 9)
        @test_throws DimensionMismatch wrap(FixedSizeArray, memoryref(m, 10), 2)
        @test wrap(FixedSizeArray, memoryref(m, 10), 1) == [10]
    end
    @testset "offset ≠ 1" begin
        for T in (Int, String)
            v = T === String ? string.(1:10) : collect(T, 1:10)
            w = copy(v)
            popfirst!(v); popfirst!(v)
            @test Base.memoryrefoffset(v.ref) == 3
            a = wrap(FixedSizeArray, v.ref, (2, 4))
            @test a == reshape(w[3:10], 2, 4)
            @test pointer(a) == pointer(v)
            @test Base.unsafe_convert(Ptr{T}, Base.cconvert(Ptr{T}, a)) == pointer(v)
            @test pointer(a, 2) == pointer(v, 2)
            a[1] = v[end]
            @test v[1] == v[end]
            @test copy(a) == a
            @test Base.memoryrefoffset(parent(copy(a))) == 1
            @test length(parent(copy(a)).mem) == 8
            @test Base.memoryrefoffset(parent(similar(a))) == 1
            @test collect(a) == reshape(v, 2, 4)
            @test Array(a) == reshape(v, 2, 4)
            @test a[:] == v
            @test vec(a) == v
            @test reshape(a, 8) == v
            @test parent(reshape(a, 4, 2)) === parent(a)
        end
        # `rand!` writes only the elements of the array
        v = zeros(10)
        popfirst!(v)
        a = wrap(FixedSizeArray, v.ref, 4)
        Random.rand!(a)
        @test all(!iszero, v[1:4])
        @test all(iszero, v[5:end])
    end
    @testset "empty" begin
        for T in (Int, String)
            v = T[]
            for size in ((0,), (0, 3), (3, 0))
                a = wrap(FixedSizeArray, v.ref, size)
                @test isempty(a)
                @test Base.size(a) === size
                @test isempty(copy(a))
                @test isempty(similar(a))
                @test collect(a) == Array{T}(undef, size)
                @test_throws BoundsErrorLight a[1]
            end
            @test wrap(MemRefFSA, v) isa FixedSizeVector{T,MemoryRef{T}}
            @test isempty(FixedSizeVector{T,MemoryRef{T}}(undef, 0))
        end
    end
    @testset "construction, `similar`, broadcasting, `copy`" begin
        for T in (Int, Float64, String)
            b = FixedSizeArray{T,2,MemoryRef{T}}(undef, 2, 3)
            @test b isa FixedSizeMatrix{T,MemoryRef{T}}
            T === String && @test !isassigned(b, 1)
            c = T === String ? string.(reshape(1:6, 2, 3)) : T.(reshape(1:6, 2, 3))
            d = MemRefFSA{T}(c)
            @test d isa FixedSizeMatrix{T,MemoryRef{T}}
            @test d == c
            @test collect_as(MemRefFSA, c) isa FixedSizeMatrix{T,MemoryRef{T}}
            @test collect_as(MemRefFSA{T,2}, c) == c
            @test similar(d) isa FixedSizeMatrix{T,MemoryRef{T}}
            @test similar(d, Int8, (4,)) isa FixedSizeVector{Int8,MemoryRef{Int8}}
            e = copy(d)
            @test e isa FixedSizeMatrix{T,MemoryRef{T}}
            @test e == d
            @test !Base.mightalias(e, d)
            e[1] = e[2]
            @test d[1] == c[1]
            if T !== String
                @test (d .+ d) isa FixedSizeMatrix{T,MemoryRef{T}}
                @test d .+ d == 2 .* c
                @test map(x -> 2x, d) == 2 .* c
                # broadcasting into an array that aliases the destination
                v = collect(T, 1:6)
                a = wrap(FixedSizeArray, v.ref, 6)
                v .= a .+ view(a, 6:-1:1)
                @test v == fill(7, 6)
                a .= v .- (1:6)
                @test v == 6:-1:1
            end
            @test FixedSizeArrays.val_parameter(FixedSizeArrays.with_stripped_type_parameters(FixedSizeArrays.TypeParametersElementType(), MemoryRef{T})) === MemoryRef
        end
        @test Base.elsize(FixedSizeMatrix{Float32,MemoryRef{Float32}}) == 4
        @test strides(FixedSizeArray{Float32,3,MemoryRef{Float32}}(undef, 2, 3, 4)) == (1, 2, 6)
    end
    @testset "`copyto!`" begin
        for T in (Int, String)
            val(i) = T === String ? string(i) : T(i)
            v = val.(1:10)
            popfirst!(v)
            a = wrap(FixedSizeArray, v.ref, 9)
            # overlapping, both directions, with the same semantics as for `Vector`
            w = copy(v)
            copyto!(a, 2, a, 1, 5)
            copyto!(w, 2, w, 1, 5)
            @test v == w
            copyto!(a, 1, a, 3, 5)
            copyto!(w, 1, w, 3, 5)
            @test v == w
            b = FixedSizeArray{T,1,MemoryRef{T}}(undef, 9)
            @test copyto!(b, a) === b
            @test b == a
            c = Vector{T}(undef, 9)
            @test copyto!(c, a) === c
            @test c == a
            @test copyto!(a, val.(11:19)) === a
            @test v == val.(11:19)
            @test copyto!(a, 3, FixedSizeVector{T}(val.(1:2)), 1, 2) === a
            @test v[3:4] == val.(1:2)
            @test_throws BoundsErrorLight copyto!(a, 9, a, 1, 2)
            @test_throws BoundsError copyto!(Vector{T}(undef, 2), a)
            @test_throws ArgumentError copyto!(a, 1, a, 1, -1)
            # different element types
            if T === Int
                f = FixedSizeArray{Float64,1,MemoryRef{Float64}}(undef, 9)
                copyto!(f, a)
                @test f == a
            end
        end
    end
    @testset "non-`isbits` elements survive GC" begin
        function make()
            v = [string(i, base = 3) for i in 1:1000]
            popfirst!(v)
            wrap(FixedSizeArray, v.ref, (999,))
        end
        a = make()
        GC.gc(true)
        @test a == [string(i, base = 3) for i in 2:1000]
        for i in eachindex(a)
            a[i] = string(i)  # write barrier: `a` is old, the strings are young
        end
        GC.gc(false); GC.gc(true)
        @test a == string.(1:999)
    end
    @testset "contract for `resize!` of the wrapped `Vector`" begin
        # the array keeps referring to the memory it was created with: memory safe,
        # but it no longer aliases the `Vector` once `resize!` moves its elements
        v = collect(1:4)
        a = wrap(FixedSizeArray, v.ref, 4)
        resize!(v, 2)  # does not move: `a` still aliases `v`
        @test Base.mightalias(a, v)
        @test a[1:2] == v
        resize!(v, 10_000)  # moves
        @test !Base.mightalias(a, v)
        @test a[1:2] == v[1:2]
        a[1] = -1
        @test v[1] == 1
        @test length(parent(a).mem) == 4
        # wrap again after resizing
        b = wrap(FixedSizeArray, v.ref, (100, 100))
        @test Base.mightalias(b, v)
        b[1] = -1
        @test v[1] == -1
        # elements of a non-`isbits` `Vector` dropped by `resize!` are unassigned
        v = ["a", "b", "c"]
        a = wrap(FixedSizeArray, v.ref, 3)
        resize!(v, 1)
        @test a[1] == "a"
        if !isassigned(a, 3)
            @test_throws UndefRefError a[3]
        end
    end
    @testset "MPI and pointer conversion" begin
        v = collect(1.0:8.0)
        popfirst!(v)
        a = wrap(FixedSizeArray, v.ref, 7)
        @test Base.unsafe_convert(MPI.MPIPtr, Base.cconvert(MPI.MPIPtr, a)) == reinterpret(MPI.MPIPtr, pointer(v))
    end
    (@isdefined AtomicMemory) && @testset "atomic memory" begin
        a = wrap(FixedSizeArray, memoryref(AtomicMemory{Int}(undef, 3)), 3)
        a .= 1:3
        @test a == 1:3
        @test copy(a) == 1:3
        @test similar(a) isa FixedSizeVector{Int,AtomicMemoryRef{Int}}
        s = wrap(FixedSizeArray, memoryref(AtomicMemory{String}(undef, 2)), 2)
        @test !isassigned(s, 1)
        s[1] = "x"
        @test s[1] == "x"
        @test copy(s)[1] == "x"
        @test !isassigned(copy(s), 2)
    end
    @testset "inference and allocations of indexing" begin
        a = FixedSizeArray{Float64,2,MemoryRef{Float64}}(undef, 3, 3)
        test_inferred(getindex, Float64, (a, 1))
        test_inferred(getindex, Float64, (a, 2, 2))
        test_noalloc(setindex!, (a, 1.0, 3))
        test_noalloc(sum, (a,))
        @test Base.infer_effects(getindex, (typeof(a), Int)) == Base.infer_effects(getindex, (FixedSizeMatrix{Float64,Memory{Float64}}, Int)) skip=!build_is_production_build
    end
    @testset "static dimensions" begin
        R{T} = FixedSizeArrayDims{T, 3, Tuple{Val{2}, Val{3}, Int}, MemoryRef{T}}
        MemRefDims = FixedSizeArrayDims{T, N, D, MemoryRef{T}} where {T, N, D}
        v = collect(1.0:30.0)
        popfirst!(v)  # a reference at an offset
        a = wrap(FixedSizeArrayDims, v.ref, (Val(2), Val(3), 4))
        @test a isa R{Float64}
        @test parent(a) === v.ref
        @test size(a) == (2, 3, 4)
        @test FixedSizeArrays.static_dims(a) === (Val(2), Val(3), 4)
        @test a == reshape(v[1:24], 2, 3, 4)
        # `:` takes all elements of the memory from the reference on
        @test size(wrap(FixedSizeArrayDims, memoryref(Memory{Float64}(undef, 30), 7), (Val(2), Val(3), :))) == (2, 3, 4)
        @test_throws DimensionMismatch wrap(FixedSizeArrayDims, v.ref, (Val(2), Val(3), :))
        # fewer elements than the memory has from the reference on
        @test size(wrap(FixedSizeArrayDims, v.ref, (Val(2), Val(3), 2))) == (2, 3, 2)
        @test_throws DimensionMismatch wrap(FixedSizeArrayDims, v.ref, (Val(2), Val(3), 5))
        @test_throws DimensionMismatch wrap(R{Float64}, v.ref, (Val(3), Val(2), 4))
        # `wrap` of a `Vector` with the type asking for the reference
        w = collect(1.0:24.0)
        b = wrap(MemRefDims, w, (Val(2), Val(3), 4))
        @test b isa R{Float64}
        @test parent(b) === w.ref
        test_inferred(wrap, R{Float64}, (MemRefDims, w, (Val(2), Val(3), 4)))
        @test b[2, 3, 4] == 24
        b[1, 2, 3] = -1
        @test w[15] == -1
        @test strides(b) == (1, 2, 6)
        @test pointer(b) == pointer(w)
        # `similar` and `copy` keep the static dimensions and the backend
        @test similar(b) isa R{Float64}
        @test similar(b, Int) isa R{Int}
        c = copy(b)
        @test c isa R{Float64}
        @test c == b
        @test parent(c) !== parent(b)
        fill!(c, 3)
        @test all(==(3), c)
        c .= b .+ 1
        @test c == b .+ 1
        @test (b .+ 1) isa FixedSizeArray{Float64, 3}
        # static sizes and strides for StaticArrayInterface.jl, and the device of the memory
        @test StaticArrayInterface.static_size(b) === (StaticInt(2), StaticInt(3), 4)
        @test StaticArrayInterface.static_strides(b) === (StaticInt(1), StaticInt(2), StaticInt(6))
        @test StaticArrayInterface.device(typeof(b)) === StaticArrayInterface.device(Memory{Float64})
        # Adapt.jl: adapting to an `Array` keeps the reference
        @test adapt(Array, b) === b
        # adapting to another storage type adapts the elements (here to a `Vector{Float32}`)
        b32 = adapt(Vector{Float32}, b)
        @test b32 isa FixedSizeArrayDims{Float32, 3, Tuple{Val{2}, Val{3}, Int}, Vector{Float32}}
        @test b32 == b
        # constructors
        @test R{Float64}(undef, (Val(2), Val(3), 5)) isa R{Float64}
        @test collect_as(R{Float64}, reshape(1.0:12.0, 2, 3, 2)) isa R{Float64}
        @test collect_as(R{Float64}, reshape(1.0:12.0, 2, 3, 2)) == reshape(1.0:12.0, 2, 3, 2)
        test_noalloc(setindex!, (b, 1.0, 2, 3, 4))
        test_inferred(getindex, Float64, (b, 2, 3, 4))
    end
end
