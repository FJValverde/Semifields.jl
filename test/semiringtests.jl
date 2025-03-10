
#using Semirings#This is already reexported in Semifields.

@testset "Boolean semiring" begin
    x, y = one(BoolSemiring), zero(BoolSemiring)
    @test valtype(typeof(x)) == Bool
    @test ! val(zero(x))
    @test val(one(x))
    #@test val(inv(one(BoolSemiring)))#Does not belong here. 
    #@test val(inv(zero(BoolSemiring)))#Does not belong here.
    @test val(x)
    @test ! val(y)
    @test val(x ⊕ y)
    @test val(x ⊕ x)
    @test ! val(y ⊕ y)
    @test ! val(x ⊗ y)
    @test val(x ⊗ x)
    @test ! val(y ⊗ y)
end

@testset "Entropy semirings" begin
    for T in [Float64, Float32]
        for a in [1, 4.5]
            K = EntropySemiring{T,a}
            x, y = K(2), K(3)
            @test valtype(typeof(x)) == T
            @test val(x ⊕ y) ≈ logaddexp(a*val(x), a*val(y)) / a
            @test val(x ⊗ y) ≈ val(x) + val(y)
            @test val(zero(x)) == T(-Inf)
            @test one(x) == K(0)
            #@test inv(x) == K(-val(x))#Does not belong here.
            #@test isnan(val(inv(zero(x))))
            @test ∂sum(x ⊕ y, x) ≈ exp(a*val(x)) / exp(a*val(x ⊕ y))
            @test ∂sum(x ⊕ y, y) ≈ exp(a*val(y)) / exp(a*val(x ⊕ y))
            @test ∂rmul(x, y) == 1
            @test ∂lmul(x, y) == 1

            z, pullback = rrule(Semirings._logaddexp , a, -Inf, -Inf)
            _, _, x̄, ȳ = pullback(1)
            @test iszero(x̄) && iszero(ȳ)
        end

        for a in [-1, -4.5]
            K = EntropySemiring{T,a}
            x, y = K(2), K(3)
            @test valtype(typeof(x)) == T
            @test val(x ⊕ y) ≈ logaddexp(a*val(x), a*val(y)) / a
            @test val(x ⊗ y) ≈ val(x) + val(y)
            @test val(zero(x)) == T(Inf)
            @test one(x) == K(0)
            #@test inv(x) ≈ K(-val(x))#Does not belong here.
            #@test isnan(val(inv(zero(x))))
            @test ∂sum(x ⊕ y, x) ≈ exp(a*val(x)) / exp(a*val(x ⊕ y))
            @test ∂sum(x ⊕ y, y) ≈ exp(a*val(y)) / exp(a*val(x ⊕ y))
            @test ∂rmul(x, y) == 1
            @test ∂lmul(x, y) == 1
        end
    end
end

#=
@testset "Probability semiring" begin
    for T in [Float32, Float64]
        K = ProbSemiring{T}
        x, y = K(2), K(3)
        @test valtype(typeof(x)) == T
        @test val(x + y) ≈ val(x) + val(y)
        @test val(x * y) ≈ val(x) * val(y)
        @test val(zero(x)) == T(0)
        @test val(one(x)) == T(1)
        #@test val(inv(x)) == inv(val(x))
        @test ∂sum(x + y, x) == one(T)
        @test ∂sum(x + y, y) == one(T)
        @test ∂rmul(x, y) == val(y)
        @test ∂lmul(x, y) == val(x)
    end
end

@testset "Log Expectation semiring" begin
    for T in [Float64, Float32]
        K = LogExpectationSemiring{T}
        semiP = EntropySemiring{T,1}
        semiV = ProbSemiring{T}
        x = K(-2.3, 0.2)
        y = K(1.3, 0.1)
        @test val(val(x + y)[1]) == val(x.prob + y.prob)
        @test val(val(x + y)[2]) == val(x.value + y.value)
        @test val(val(x * y)[1]) == val(x.prob * y.prob)
        @test val(val(x * y)[2]) ≈ val((semiV(exp(x.prob.val)) * y.value) + (semiV(exp(y.prob.val)) * x.value) )
        @test val(zero(x)) == (zero(semiP), zero(semiV))
        @test val(one(x)) == (one(semiP), one(semiV))
    end
end

@testset "Prob Expectation Semiring" begin
    for T in [Float64, Float32]
            K = ProbExpectationSemiring{T}
			S = ProbSemiring{T}
			x = K(-2.3, 0.2)
			y = K(1.3, 0.1)
            @test val(val(x + y)[1]) == val(x.prob + y.prob)
			@test val(val(x + y)[2]) == val(x.value + y.value)
            @test val(val(x * y)[1]) == val(x.prob * y.prob)
			@test val(val(x * y)[2]) ≈ val((S(identity(x.prob.val)) * y.value) + (S(identity(y.prob.val)) * x.value) )
            @test val(zero(x)) == (zero(S), zero(S))
            @test val(one(x)) == (one(S), one(S))
    end
end
=#

@testset "Tropical/Arctic semiring" begin
    for T in [Float32, Float64]
        K = TropicalSemiring{T}
        x, y = K(2), K(3)
        @test valtype(typeof(x)) == T
        @test val(x ⊕ y) ≈ min(val(x), val(y))
        @test val(x ⊗ y) ≈ val(x) + val(y)
        @test val(zero(x)) == T(Inf)
        @test val(one(x)) == T(0)
        #@test val(inv(x)) == T(-val(x))
        #@test isnan(val(inv(zero(x))))
        #@test ∂sum(x ⊕ y, x) ≈ T(val(x) ≤ val(y))#Don't like the implicit conversion.
        #@test ∂sum(x ⊕ y, y) ≈ T(val(y) ≤ val(x))#Don't like the implicit conversion.
        @test ∂rmul(x, y) == 1
        @test ∂lmul(x, y) == 1

        K = ArcticSemiring{T}
        x, y = K(2), K(3)
        @test valtype(typeof(x)) == T
        @test val(x ⊕ y) ≈ max(val(x), val(y))
        @test val(x ⊗ y) ≈ val(x) + val(y)
        @test val(zero(x)) == T(-Inf)
        @test val(one(x)) == T(0)
        #@test val(inv(x)) == T(-val(x))
        #@test isnan(val(inv(zero(x))))
        #@test ∂sum(x ⊕ y, x) ≈ T(val(x) ≥ val(y))#Don't like the implicit conversion.
        #@test ∂sum(x ⊕ y, y) ≈ T(val(y) ≥ val(x))#Don't like the implicit conversion.
        @test ∂rmul(x, y) == 1
        @test ∂lmul(x, y) == 1
    end
end

@testset "conversion" begin
    for S in [TropicalSemiring{Float32}, EntropySemiring{Float64, -1}]
        @test all(S[1, 2, 3] .== [S(1), S(2), S(3)])
    end

    @test convert(TropicalSemiring{Float32}, one(EntropySemiring{Float32,1})) isa TropicalSemiring{Float32}
    @test convert(TropicalSemiring{Float32}, one(ProbSemiring{Float64})) isa TropicalSemiring{Float32} skip = true
    @test convert(BoolSemiring, one(EntropySemiring{Float32,1})) isa BoolSemiring
    @test !isone(convert(BoolSemiring, one(EntropySemiring{Float32,1})))
    @test convert(BoolSemiring, zero(EntropySemiring{Float32,1})) isa BoolSemiring skip=true
    @test iszero(convert(BoolSemiring, zero(EntropySemiring{Float32,1}))) skip=true
end

#=FVA: there is no need for this. TO be solved on semimodules. 
@testset "integer multiplication" begin
    for S in [TropicalSemiring{Float32}, EntropySemiring{Float64, 2.1}]
        x = S(2.3)

        @test val(3 ⊗ x) ≈ val(x ⊗ x ⊗ x)
        @test val(x ⊗ 4) ≈ val(x ⊗ x ⊗ x ⊗ x)
        @test iszero(0 ⊗ x)
    end
end
=#

@testset "morphisms" begin
	for T in [Float64, Float32]
        x = EntropySemiring{T,1}(2.3)
        y = EnergySemiring{T}(exp(2.3))
        @test val(x) ≈ log(val(y))
        @test val(y) ≈ exp(val(x))
	end
end


