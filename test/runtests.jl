using Test
using FixedSizeArrays
using OffsetArrays: OffsetArray
import Aqua

# Check if the compilation options allow maximum performance.
const build_is_production_build_env_name = "BUILD_IS_PRODUCTION_BUILD"
const build_is_production_build = let v = get(ENV, build_is_production_build_env_name, "true")
    if v ∉ ("false", "true")
        error("unknown value for environment variable $build_is_production_build_env_name: $v")
    end
    if v == "true"
        true
    else
        false
    end
end::Bool

const checked_dims = FixedSizeArrays.checked_dims

# helpers for testing for allocation or suboptimal inference

allocated(f::F,      args::Tuple) where {F} = @allocated f(args...)
allocated(::Type{T}, args::Tuple) where {T} = @allocated T(args...)
test_noalloc(f::F,      args::Tuple) where {F} = build_is_production_build && @test iszero(allocated(f, args))
test_noalloc(::Type{T}, args::Tuple) where {T} = build_is_production_build && @test iszero(allocated(T, args))
function test_inferred(f::F,      ::Type{R}, args::Tuple) where {F,R}
    @test isconcretetype(R)  # meta: test-strictness test
    @test (@inferred f(args...)) isa R
end
function test_inferred(::Type{T}, ::Type{R}, args::Tuple) where {T,R}
    @test isconcretetype(R)  # meta: test-strictness test
    @test (@inferred T(args...)) isa R
end
test_inferred(::Type{T}, args::Tuple) where {T} = test_inferred(T, T, args)
function test_inferred_noalloc(f::F,      ::Type{R}, args::Tuple) where {F,R}
    test_noalloc(f, args)
    test_inferred(f, R, args)
end
function test_inferred_noalloc(::Type{T}, ::Type{R}, args::Tuple) where {T,R}
    test_noalloc(T, args)
    test_inferred(T, R, args)
end
function test_inferred_noalloc(::Type{T},            args::Tuple) where {T  }
    test_noalloc(T, args)
    test_inferred(T, args)
end

# helpers for constructing the type constructors

function fsa(vec_type::Type{<:DenseVector})
    FixedSizeArray{T,N,vec_type{T}} where {T,N}
end
function fsm(vec_type::Type{<:DenseVector})
    fsa(vec_type){T,2} where {T}
end
function fsv(vec_type::Type{<:DenseVector})
    fsa(vec_type){T,1} where {T}
end

function test_we_do_not_own_the_call(func, arg_types)
    function f(method)
        FixedSizeArrays != parentmodule(method)
    end
    @test all(f, methods(func, arg_types))
end

@testset "FixedSizeArrays" begin
    @testset "meta" begin
        @test isempty(detect_ambiguities(Main))  # test for ambiguities in this file
    end

    @testset "Aqua.jl" begin
        Aqua.test_all(FixedSizeArrays)
    end

    @testset "type piracy" begin
        @testset "issue #77: type piracy of `similar`" begin
            test_we_do_not_own_the_call(similar, Tuple{Broadcast.Broadcasted{Broadcast.ArrayStyle{Union{}}}, Type{Int}})
        end
    end

    @testset "type aliases" begin
        @test FixedSizeArrayDefault <: FixedSizeArray
        @test FixedSizeVectorDefault <: FixedSizeVector
        @test FixedSizeMatrixDefault <: FixedSizeMatrix
        for et ∈ (Float32, Int, String)
            @test FixedSizeArrayDefault{et} <: FixedSizeArray{et}
            @test FixedSizeVectorDefault{et} <: FixedSizeVector{et}
            @test FixedSizeMatrixDefault{et} <: FixedSizeMatrix{et}
            @test FixedSizeArrayDefault{et, 0} <: FixedSizeArray{et, 0}
            for t ∈ (
                FixedSizeVectorDefault{et}, FixedSizeMatrixDefault{et}, FixedSizeArrayDefault{et, 0},
            )
                @test isconcretetype(t)
            end
        end
    end

    @testset "safe computation of length from dimensions size" begin
        @test isone(checked_dims(()))
        for n ∈ 0:30
            t = Tuple(1:n)
            if 20 < n
                @test_throws ArgumentError checked_dims(t)
            else
                @test factorial(n) == prod(t) == checked_dims(t)
                test_inferred_noalloc(checked_dims, Int, (t,))
            end
        end
    end

    struct Iter{E,N,I<:Integer}
        size::NTuple{N,I}
        length::I
        val::E
    end
    function Base.iterate(i::Iter)
        l = i.length
        iterate(i, max(zero(l), l))
    end
    function Base.iterate(i::Iter, state::Int)
        if iszero(state)
            nothing
        else
            (i.val, state - 1)
        end
    end
    Base.IteratorSize(::Type{<:Iter{<:Any,N}}) where {N} = Base.HasShape{N}()
    Base.length(i::Iter) = i.length
    Base.size(i::Iter) = i.size
    Base.eltype(::Type{<:Iter{E}}) where {E} = E

    @testset "default underlying storage type" begin
        default = FixedSizeArrays.default_underlying_storage_type
        @test default === (@isdefined(Memory) ? Memory : Vector)
        return_type = FixedSizeVector{Int,default{Int}}
        test_inferred(FixedSizeArray{Int}, return_type, (undef, 3))
        test_inferred(FixedSizeArray{Int}, return_type, (undef, (3,)))
        test_inferred(FixedSizeVector{Int}, return_type, (undef, 3))
        test_inferred(FixedSizeVector{Int}, return_type, (undef, (3,)))
        iter = Iterators.filter(iseven, 3:7)
        @test collect_as(FixedSizeArray, iter) isa return_type
        test_inferred(collect_as, return_type, (FixedSizeArray{Int}, iter))
        test_inferred(collect_as, return_type, (FixedSizeVector{Int}, iter))
        arr = ([1, 2, 3],)
        test_inferred(FixedSizeArray, return_type, arr)
        test_inferred(FixedSizeVector, return_type, arr)
        test_inferred(FixedSizeArray{Int}, return_type, arr)
        test_inferred(FixedSizeVector{Int}, return_type, arr)
    end

    for storage_type ∈ (((@isdefined Memory) ? (Memory,) : ())..., Vector)
        FSV = fsv(storage_type)
        FSM = fsm(storage_type)
        FSA = fsa(storage_type)

        @testset "`propertynames`" begin
            for dim_count ∈ 0:4
                for T ∈ (Int, Float32)
                    siz = ntuple(Returns(1), dim_count)
                    arr = FSA{T}(undef, siz)
                    test_inferred(propertynames, Tuple{}, (arr,))
                    test_inferred(propertynames, Tuple{}, (arr, false))
                end
            end
        end

        @testset "Constructors" begin
            for dim_count ∈ 0:4
                siz = ntuple(Returns(2), dim_count)
                T = FSA{Float64}
                R = FSA{Float64,dim_count}
                for sz ∈ (siz, map(BigInt, siz))
                    for args ∈ ((undef, sz), (undef, sz...))
                        test_inferred(   R, args)
                        test_inferred(T, R, args)
                    end
                end
            end
            for offset ∈ (-1, 0, 2, 3)
                ax = offset:(offset + 1)
                oa = OffsetArray([10, 20], ax)
                @test_throws DimensionMismatch FSA{Int}(oa)
                @test_throws DimensionMismatch FSV{Int}(oa)
            end
            @testset "taken from Julia's test/core.jl" begin
                # inspired by:
                # https://github.com/JuliaLang/julia/blob/83929ad883c97fa4376618c4559b74f6ada7a3ce/test/core.jl#L7211-L7228
                b = prevpow(2, typemax(Int))
                test_inferred(FSA{Int}, FSA{Int,3}, (undef, 0, b, b))
                test_inferred(FSA{Int}, FSA{Int,3}, (undef, b, b, 0))
                test_inferred(FSA{Int}, FSA{Int,5}, (undef, b, b, 0, b, b))
                @test_throws ArgumentError FSA{Int}(undef, b, b)
                @test_throws ArgumentError FSA{Int}(undef, 1, b, b)
                @test_throws ArgumentError FSA{Int}(undef, 0, -10)
                @test_throws ArgumentError FSA{Int}(undef, -10, 0)
                @test_throws ArgumentError FSA{Int}(undef, -1, -1)
                @test_throws ArgumentError FSA{Int}(undef, 0, -4, -4)
                @test_throws ArgumentError FSA{Int}(undef, -4, 1, 0)
                @test_throws ArgumentError FSA{Int}(undef, -4, -4, 1)
            end
            @test_throws ArgumentError FSA{Float64,1}(undef, -1)
            @test_throws ArgumentError FSA{Float64,1}(undef, (-1,))
            @test_throws ArgumentError FSA{Float64,2}(undef, -1, -1)
            @test_throws ArgumentError FSA{Float64,2}(undef, typemax(Int), 0)
            @test_throws ArgumentError FSA{Float64,2}(undef, typemax(Int), typemax(Int))
            @test_throws ArgumentError FSA{Float64,3}(undef, typemax(Int)-1, typemax(Int)-1, 2)
            @test_throws ArgumentError FSA{Float64,4}(undef, typemax(Int)-1, typemax(Int)-1, 2, 4)
            @testset "negative dimension size" begin
                for n ∈ 2:4
                    funs = (
                        Returns(-1),
                        (i -> (i == 1) ?  1 : -1),
                        (i -> (i == 1) ? -1 :  1),
                        (i -> (i == n) ?  1 : -1),
                        (i -> (i == n) ? -1 :  1),
                    )
                    fun = f -> ntuple(f, n)
                    sizes = map(fun, funs)
                    for siz ∈ sizes
                        @test_throws ArgumentError FSA{Float64,n}(undef, siz)
                        @test_throws ArgumentError FSA{Float64,n}(undef, siz...)
                    end
                end
            end
        end

        @testset "setindex!" begin
            v = FSV{Float64}(undef, 2)
            v[1] = 1
            v[2] = 2
            @test v[1] == 1
            @test v[2] == 2
            @test_throws BoundsError v[3] = 3
        end

        @testset "FixedSizeVector" begin
            v = FSV{Float64}(undef, 3)
            @test length(v) == 3
            test_inferred_noalloc(((a, b) -> a .= b), FSV{Float64}, (v, 1:3))
            @test v == 1:3
            test_inferred(similar, FSV{Float64}, (v,))
            test_inferred(similar, FSV{Int}, (v, Int))
            test_inferred(copy, FSV{Float64}, (v,))
            test_inferred(zero, FSV{Float64}, (v,))
            test_inferred(similar, FSV{Int}, (FSV{Int}, (2,)))
            test_inferred(similar, FSV{Int}, (FSA{Int}, (2,)))
            example_abstract_vectors = (7:9, [7, 8, 9], OffsetArray([7, 8, 9], 1:3))
            test_convert = function(T, a)
                test_inferred(convert, FSV{Int}, (T, a))
                c = convert(T, a)
                test_inferred(convert, FSV{Int}, (T, c))
                @test c == a
                @test c === convert(T, c)
            end
            for T ∈ (FSA{Int}, FSV{Int})
                for S ∈ (Int, Float64)
                    for a ∈ map((c -> map(S, c)), example_abstract_vectors)
                        test_convert(T, a)
                    end
                end
            end
        end

        @testset "FixedSizeMatrix" begin
            m = FSM{Float64}(undef, 3, 3)
            m_ref = reshape(1:9, size(m))
            @test length(m) == 9
            test_inferred_noalloc(((a, b) -> a .= b), FSM{Float64}, (m, m_ref))
            @test m == m_ref
            test_inferred(similar, FSM{Float64}, (m,))
            test_inferred(similar, FSM{Int}, (m, Int))
            test_inferred(copy, FSM{Float64}, (m,))
            test_inferred(zero, FSM{Float64}, (m,))
            test_inferred(similar, FSM{Int}, (FSM{Int}, (2, 3)))
            test_inferred(similar, FSM{Int}, (FSA{Int}, (2, 3)))
            example_abstract_matrices = (
                reshape(1:9, (3, 3)),
                OffsetArray(reshape(1:9, (3, 3)), 1:3, 1:3),
            )
            test_convert = function(T, a)
                test_inferred(convert, FSM{Int}, (T, a))
                c = convert(T, a)
                test_inferred(convert, FSM{Int}, (T, c))
                @test c == a
                @test c === convert(T, c)
            end
            for T ∈ (FSA{Int}, FSM{Int})
                for S ∈ (Int, Float64)
                    for a ∈ map((c -> map(S, c)), example_abstract_matrices)
                        test_convert(T, a)
                    end
                end
            end
        end

        @testset "`map`" begin
            for dim_count ∈ 0:3
                size = ntuple(Returns(3), dim_count)
                a = FSA{Int, dim_count}(undef, size)
                for v ∈ (3, 3.1, nothing)
                    args = (Returns(v), a)
                    @test all(==(v), map(args...))
                    test_inferred(map, FSA{typeof(v), dim_count}, args)
                end
            end
        end

        @testset "broadcasting" begin
            TVI = FSV{Int}
            TVF = FSV{Float64}
            v3 = TVI(1:3)
            vf3 = TVF(1:3)
            bp = (x, y) -> x .+ y
            bm = (x, y) -> x .* y
            test_inferred(+,  TVI, (v3, v3))
            test_inferred(bp, TVI, (v3, v3))
            test_inferred(bm, TVI, (v3, v3))
            test_inferred(+,  TVF, (v3, vf3))
            test_inferred(bp, TVF, (v3, vf3))
            test_inferred(bm, TVF, (v3, vf3))
            test_inferred(bp, TVI, (v3,  3))
            test_inferred(bm, TVI, (v3,  3))
            test_inferred(bp, TVF, (v3, .3))
            test_inferred(bm, TVF, (v3, .3))
            @testset "matrices" begin
                TMI = FSM{Int}
                TMF = FSM{Float64}
                m33 = TMI(undef, 3, 3)
                m13 = TMI(undef, 1, 3)
                m31 = TMI(undef, 3, 1)
                test_inferred(+,  TMI, (m33, m33))
                test_inferred(bp, TMI, (m33, m33))
                test_inferred(bm, TMI, (m33, m33))
                test_inferred(bp, TMI, (m33, m13))
                test_inferred(bm, TMI, (m33, m13))
                test_inferred(bp, TMI, (m33, m31))
                test_inferred(bm, TMI, (m33, m31))
                test_inferred(bp, TMI, (m33, v3 ))
                test_inferred(bm, TMI, (m33, v3 ))
                test_inferred(bp, TMF, (m33, vf3))
                test_inferred(bm, TMF, (m33, vf3))
                test_inferred(bp, TMI, (m33,   3))
                test_inferred(bm, TMI, (m33,   3))
                test_inferred(bp, TMF, (m33,  .3))
                test_inferred(bm, TMF, (m33,  .3))
            end
        end

        @testset "`isassigned`" begin
            for dim_count ∈ 0:3
                for elem_type ∈ (Int, FSM{Nothing})
                    size = ntuple(Returns(3), dim_count)
                    a = FSA{elem_type, dim_count}(undef, size)
                    for i_style ∈ (IndexLinear(), IndexCartesian())
                        for i ∈ eachindex(i_style, a)
                            @test isassigned(a, i) == isbitstype(elem_type)
                            test_inferred_noalloc(isassigned, Bool, (a, i))
                        end
                    end
                end
            end
            @testset "some assigned" begin
                a = FSM{FSM{Nothing}}(undef, 3, 3)
                assigned_inds = ((1, 2), (2, 3), (3, 3))
                for ij ∈ assigned_inds
                    a[ij...] = FSM{Nothing}(undef, 1, 1)
                end
                for ij ∈ Iterators.product(1:3, 1:3)
                    @test isassigned(a, ij...) == (ij ∈ assigned_inds)
                    test_inferred_noalloc(isassigned, Bool, (a, ij...))
                end
            end
        end

        @testset "`copyto!`" begin
            for (D, S) ∈ (
                (Vector, FSV),
                (FSV, Vector),
                (FSV, FSV),
                (
                    if @isdefined Memory
                        ((FSV, Memory), (Memory, FSV))
                    else
                        ()
                    end
                )...,
            )
                for E ∈ (Float64, Int)
                    s = S{E}(undef, 5)
                    s .= 1:5
                    d = D{Float64}(undef, 5)
                    test_inferred_noalloc(copyto!, D{Float64}, (d, s))
                    @test copyto!(d, s) == 1:5
                    args = (d, 1, s, 1, length(s))
                    test_inferred_noalloc(copyto!, D{Float64}, args)
                    @test copyto!(args...) == 1:5
                    @test_throws ArgumentError copyto!(d, 1, s, 1, -1)
                end
            end
        end

        @testset "LinearAlgebra" begin
            @testset "$T" for T in (Float16, Float32, Float64)
                v = randn(T, 8)
                v_fixed = FSV{T}(v)
                test_inferred(+, FSV{T}, (v_fixed, v_fixed))
                v_sum = v_fixed + v_fixed
                @test v_sum ≈ v + v
                test_inferred(((x, y) -> x * y'), FSM{T}, (v_fixed, v_fixed))
                v_mul = v_fixed * v_fixed'
                @test v_mul ≈ v * v'

                M = randn(T, 8, 8)
                M_fixed = FSM{T}(M)
                test_inferred(+, FSM{T}, (M_fixed, M_fixed))
                M_sum = M_fixed + M_fixed
                @test M_sum ≈ M + M
                test_inferred(*, FSM{T}, (M_fixed, M_fixed))
                M_mul = M_fixed * M_fixed
                @test M_mul ≈ M * M
            end
        end

        @testset "`pointer`" begin
            for elem_type ∈ (Int8, Int16, Int32, Int64)
                for dim_count ∈ 0:4
                    siz = ntuple(Returns(2), dim_count)
                    arr = FSA{elem_type, dim_count}(undef, siz)
                    @test pointer(arr) === pointer(arr, 1) === pointer(arr, 2) - sizeof(elem_type)
                    test_inferred_noalloc(pointer, Ptr{elem_type}, (arr,))
                    test_inferred_noalloc(pointer, Ptr{elem_type}, (arr, 1))
                end
            end
        end

        @testset "`reshape`" begin
            length_to_shapes = Dict(
                (0 => ((0,), (0, 0), (0, 1), (1, 0), (1, 0, 0), (0, 0, 1))),
                (1 => ((), (1,), (1, 1), (1, 1, 1))),
                (2 => ((2,), (1, 2), (2, 1), (1, 2, 1))),
                (3 => ((3,), (1, 3), (3, 1), (1, 3, 1))),
                (4 => ((4,), (1, 4), (4, 1), (2, 2), (1, 2, 2), (2, 1, 2))),
                (6 => ((6,), (1, 6), (6, 1), (2, 3), (3, 2), (1, 3, 2), (2, 1, 3))),
            )
            for elem_type ∈ (Int, Number, Union{Nothing,Int})
                for len ∈ keys(length_to_shapes)
                    shapes = length_to_shapes[len]
                    for shape1 ∈ shapes
                        a = FSA{elem_type,length(shape1)}(undef, shape1)
                        @test_throws DimensionMismatch reshape(a, length(a)+1)
                        @test_throws DimensionMismatch reshape(a, length(a)+1, 1)
                        @test_throws DimensionMismatch reshape(a, 1, length(a)+1)
                        for shape2 ∈ shapes
                            @test prod(shape1) === prod(shape2) === len  # meta
                            T = FSA{elem_type,length(shape2)}
                            test_inferred_noalloc(reshape, T, (a, shape2))
                            test_inferred_noalloc(reshape, T, (a, shape2...))
                            b = reshape(a, shape2)
                            @test size(b) === shape2
                            @test a.mem === b.mem
                            @test a === reshape(b, shape1)
                        end
                    end
                end
            end
        end

        @testset "`collect_as`" begin
            @testset "difficult requested return type" begin
                T = FixedSizeVectorDefault{T} where {T <: Float32}
                iter = 3:7
                # either return a value of the requested type or throw
                returns = try
                    collect_as(T, iter)
                    true
                catch e
                    (e isa TypeError) || rethrow()
                    false
                end
                returns && @test collect_as(T, iter) isa T
            end
            @testset "empty iterator with inexact `eltype`" begin
                iterator = Iterators.map((x -> x + 0.3), [])
                @test collect_as(FSV, iterator) isa FSV{Union{}}
                @test collect_as(FSV{Union{}}, iterator) isa FSV{Union{}}
                @test collect_as(FSV{Float32}, iterator) isa FSV{Float32}
                @test collect_as(FSV{Any}, iterator) isa FSV{Any}
            end
            @testset "`Union{}`" begin
                @test_throws Exception collect_as(Union{}, ())
                @test_throws Exception collect_as(FixedSizeVector{Union{}, Union{}}, ())
                @test_throws Exception collect_as(FixedSizeVector{<:Any, Union{}}, ())
            end
            for T ∈ (FSA{Int}, FSV{Int})
                for iterator ∈ (Iterators.repeated(7), Iterators.cycle(7))
                    @test_throws ArgumentError collect_as(T, iterator)
                end
            end
            @test_throws ArgumentError collect_as(FSA{Int, -1}, 7:8)
            @test_throws TypeError collect_as(FSA{Int, 3.1}, 7:8)
            for T ∈ (FSA{3}, FSV{3})
                iterator = (7:8, (7, 8))
                @test_throws TypeError collect_as(T, iterator)
            end
            @testset "buggy iterator with mismatched `size` and `length" begin
                for iterator ∈ (Iter((), 0, 7), Iter((3, 2), 5, 7))
                    E = eltype(iterator)
                    dim_count = length(size(iterator))
                    for T ∈ (FSA{E}, FSA{E,dim_count})
                        @test_throws DimensionMismatch collect_as(T, iterator)
                    end
                end
            end
            iterators = (
                (), (7,), (7, 8), 7, (7 => 8), Ref(7), fill(7),
                (i for i ∈ 1:3), ((i + 100*j) for i ∈ 1:3, j ∈ 1:2), Iterators.repeated(7, 2),
                (i for i ∈ 7:9 if i==8), 7:8, 8:7, map(BigInt, 7:8), Int[], [7], [7 8],
                Iter((), 1, 7), Iter((3,), 3, 7), Iter((3, 2), 6, 7),
                (i for i ∈ (false, 0x1, 2) if Bool(2 - 1)),
                (i + false for i ∈ (false, 0x1, 2) if Bool(2 - 2)),
                Iterators.filter(<(1), Int[]),
            )
            abstract_array_params(::AbstractArray{T,N}) where {T,N} = (T, N)
            @testset "iterator: $iterator" for iterator ∈ iterators
                a = collect(iterator)
                (E, dim_count) = abstract_array_params(a)
                af = collect(Float64, iterator)
                @test abstract_array_params(af) == (Float64, dim_count)  # meta
                @test_throws DimensionMismatch collect_as(FSA{E,dim_count+1}, iterator)
                for T ∈ (FSA{E}, FSA{E,dim_count})
                    test_inferred(collect_as, FSA{E,dim_count}, (T, iterator))
                    fsa = collect_as(T, iterator)
                    @test a == fsa
                    @test first(abstract_array_params(fsa)) <: E
                end
                for T ∈ (FSA, FSA{<:Any,dim_count})
                    @test collect_as(T, iterator) isa FSA{E,dim_count}
                    fsa = collect_as(T, iterator)
                    @test a == fsa
                    @test first(abstract_array_params(fsa)) <: E
                end
                for T ∈ (FSA{Float64}, FSA{Float64,dim_count})
                    test_inferred(collect_as, FSA{Float64,dim_count}, (T, iterator))
                    fsa = collect_as(T, iterator)
                    @test af == fsa
                    @test first(abstract_array_params(fsa)) <: Float64
                end
                for T ∈ (FixedSizeArray{<:Any,dim_count},)
                    @test collect_as(T, iterator) isa FixedSizeArray{E,dim_count}
                    fsa = collect_as(T, iterator)
                    @test a == fsa
                    @test first(abstract_array_params(fsa)) <: E
                end
            end
        end
    end
end
