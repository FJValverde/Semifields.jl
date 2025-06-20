#using Semifields
#using Test
#using ChainRulesCore
import LogExpFunctions: logaddexp
@testset "Boolean semifield" begin
    x, y = one(BoolSemifield), zero(BoolSemifield)
    for v in [true, false]
        @test val(BoolSemifield(v)) == v
    end
    @test valtype(typeof(x)) == Bool
    @test ! val(zero(x))
    @test val(one(x))
    @test val(x)
    @test ! val(y)
    @test val(x ⊕ y)
    @test val(x ⊕ x)
    @test ! val(y ⊕ y)
    @test ! val(x ⊗ y)
    @test val(x ⊗ x)
    @test ! val(y ⊗ y)
end

#The following behaviour lies at the centre of every complete semifield.
@testset "zero/one elements in Entropy Semifields" begin
    for T in [Float64, Float32]
        for S in [#ProbSemifield{T},
                  EntropySemifield{T, -0.001},#Shannon semifield, approx
                  EntropySemifield{T, 0.001}, #Dual Shannon semifield, approx
                  EntropySemifield{T, -1},#harmonic semifield
                  EntropySemifield{T, 1}, #standard semifield
                  TropicalSemifield{T},
                  ArcticSemifield{T},
                  MaxPlusSemifield{T},
                  MinPlusSemifield{T}
                  ]
            for v in [-Inf, -1.0, 0.0, 1.0, Inf]
                @test val(S(v)) == v
            end
            #@test inv(zero(S)) == top(S)
            #@test inv(top(S)) == zero(S)
            #@test inv(one(S)) ≈ one(S)#returns -0.0 != 0.0
            for y in [zero(S), one(S)]#, top(S)]
                @test zero(S) ⊗ y == zero(S)
                @test y ⊗ zero(S) == zero(S)
                @test zero(S) ⊕ y == y
                @test y ⊕ zero(S) == y
                @test one(S) ⊗ y == y
                @test y ⊗ one(S) == y
                #@test top(S)  ⊕ y == top(S)
                #@test y  ⊕ top(S) == top(S)
            end
        end
    end
end

@testset "Entropy semifield operations" begin
    for T in [Float64, Float32]
        for a in [1.0, 4.5]
            K = EntropySemifield{T,a}
            x, y = K(2), K(3)
            @test valtype(typeof(x)) == T
            @test val(x ⊕ y) ≈ logaddexp(a*val(x), a*val(y)) / a
            @test val(x ⊗ y) ≈ val(x) + val(y)
            @test val(zero(x)) == T(-Inf)
            @test val(one(x)) == T(0)
            @test ∂sum(x ⊕ y, x) ≈ exp(a*val(x)) / exp(a*val(x ⊕ y))
            @test ∂sum(x ⊕ y, y) ≈ exp(a*val(y)) / exp(a*val(x ⊕ y))
            @test ∂rmul(x, y) == 1
            @test ∂lmul(x, y) == 1

            z, pullback = rrule(Semifields._logaddexp , a, -Inf, -Inf)
            _, _, x̄, ȳ = pullback(1)
            @test iszero(x̄) && iszero(ȳ)
        end

        for a in [-1.0, -4.5]
            K = EntropySemifield{T,a}
            x, y = K(2), K(3)
            @test valtype(typeof(x)) == T
            @test val(x ⊕ y) ≈ logaddexp(a*val(x), a*val(y)) / a
            @test val(x ⊗ y) ≈ val(x) + val(y)
            @test val(zero(x)) == T(Inf)
            @test val(one(x)) == T(0)
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
        K = ProbSemifield{T}
        x, y = K(2), K(3)
        @test valtype(typeof(x)) == T
        @test val(x ⊕ y) ≈ val(x) + val(y)
        @test val(x ⊗ y) ≈ val(x) * val(y)
        @test val(zero(x)) == T(0)
        @test val(one(x)) == T(1)
        @test inv(zero(K)) == top(K)
        @test inv(top(K)) == zero(K)
        @test ∂sum(x ⊕ y, x) ≈ 1
        @test ∂sum(x ⊕ y, y) ≈ 1
        @test ∂rmul(x, y) == val(y)
        @test ∂lmul(x, y) == val(x)
    end
end
=#
@testset "Tropical/Arctic semiring" begin
    for T in [Float32, Float64]
        for K in [TropicalSemifield{T}, MinPlusSemifield{T}]
        x, y = K(2), K(3)
        @test valtype(typeof(x)) == T
        @test val(x ⊕ y) ≈ min(val(x), val(y))
        @test val(x ⊗ y) ≈ val(x) + val(y)
        @test val(zero(x)) == T(Inf)
        @test val(one(x)) == T(0)
        @test ∂sum(x ⊕ y, x) ≈ T(val(x) ≤ val(y))
        @test ∂sum(x ⊕ y, y) ≈ T(val(y) ≤ val(x))
        @test ∂rmul(x, y) == 1
        @test ∂lmul(x, y) == 1
        end

        for K in [ArcticSemifield{T}, MaxPlusSemifield{T}]
        x, y = K(2), K(3)
        @test valtype(typeof(x)) == T
        @test val(x ⊕ y) ≈ max(val(x), val(y))
        @test val(x ⊗ y) ≈ val(x) + val(y)
        @test val(zero(x)) == T(-Inf)
        @test val(one(x)) == T(0)
        @test ∂sum(x ⊕ y, x) ≈ T(val(x) ≥ val(y))
        @test ∂sum(x ⊕ y, y) ≈ T(val(y) ≥ val(x))
        @test ∂rmul(x, y) == 1
        @test ∂lmul(x, y) == 1
        end
    end
end

@testset "conversion" begin
    for T in [Float32, Float64]
        for S in [TropicalSemifield{T}, MinPlusSemifield{T},
                    ArcticSemifield{T}, MaxPlusSemifield{T},
                    EntropySemifield{T, -1},
                    EntropySemifield{T, 1},
                    EntropySemifield{T, -0.001},
                    EntropySemifield{T, 0.001}]
            @test all(S[1, 2, 3] .== [S(1), S(2), S(3)])
        end
    end
end

#= FVA THese have no sense for the time being. 
@testset "integer multiplication" begin
    for S in [TropicalSemifield{Float32}, EntropySemifield{Float64, 2.1}]
        x = S(2.3)

        @test val(3 * x) ≈ val(x ⊕ x ⊕ x)
        @test val(x * 4) ≈ val(x ⊕ x ⊕ x ⊕ x)
        @test iszero(0 * x)
    end
end
=#