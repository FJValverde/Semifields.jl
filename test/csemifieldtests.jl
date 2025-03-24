#using Semifields
#using Test
#using ChainRulesCore
#import LogExpFunctions: logaddexp
using ..CSemifields

@testset "Complete Semifields: unitary constructors" begin
    for S in [TernaryCSemifield]
        @test ⊤(S) == top(S) && ⊥(S) == zero(S) && ⌶(S) == one(S)
        @test inv(⌶(S)) == ⌶(S)
        @test inv(⊤(S)) == ⊥(S)
        @test inv(⊥(S)) == ⊤(S)
    end

    for T in [Float32, Float64]
        for τ in [-Inf, -1.0, -1e-10, 1e-10, 1.0, Inf]
            for S in [EntropyCSemifield{T,τ}]
                @test ⊤(S) == top(S) && ⊥(S) == zero(S) && ⌶(S) == one(S)
                @test inv(⌶(S)) == ⌶(S)
                @test inv(⊤(S)) == ⊥(S)
                @test inv(⊥(S)) == ⊤(S)
            end
        end
        for S in [TropicalCSemifield{T}, MinPlusCSemifield{T}, 
                  ArcticCSemifield{T}, MaxPlusCSemifield{T}]
            @test ⊤(S) == top(S) && ⊥(S) == zero(S) && ⌶(S) == one(S)
            @test inv(⌶(S)) == ⌶(S)
            @test inv(⊤(S)) == ⊥(S)
            @test inv(⊥(S)) == ⊤(S)
        end
    end
end
import ..Semirings: _logaddexp

@testset "Complete Semifields: basic operations" begin
    for T in [Float32, Float64]
        for τ in [-10, -1.0, -1e-10, 1e-10, 1.0, 10]
            for S in [EntropyCSemifield{T,τ}]
                x = S(1.0)
                y = S(2.0)
                z = S(3.0)
                @test x ⊕ y == S( _logaddexp(τ, val(x), val(y)) )
                @test x ⊗ y == S(3.0)
                @test x ⊗ ⌶(S) == x
                @test x ⊕ ⊤(S) == ⊤(S)
                @test x ⊗ ⊥(S) == ⊥(S)
                @test ⊤(S) ⊕ ⊤(S) == ⊤(S)#idempotence of addition
                @test ⊤(S) ⊗ ⊤(S) == ⊤(S)#idempotence of multiplication
                @test ⊤(S) ⊕ ⊥(S) == ⊤(S)#null element for addition
                @test ⊤(S) ⊗ ⊥(S) == ⊥(S)#absorptive element for multiplication
                @test ⊥(S) ⊕ ⊤(S) == ⊤(S)# null element for addition, commutativity
                @test ⊥(S) ⊗ ⊤(S) == ⊥(S)#absorptive element for multiplication, commutativity
                @test ⊥(S) ⊕ ⌶(S) == ⌶(S)#identity element for addition
                @test ⊥(S) ⊗ ⌶(S) == ⊥(S)#identity element for multiplication
                @test ⌶(S) ⊕ ⊤(S) == ⊤(S)#absorptive element for addition
                @test ⌶(S) ⊗ ⊤(S) == ⊤(S)#identity element for multiplication
                @test ⌶(S) ⊕ ⊥(S) == ⌶(S)#identity element for addition
                @test ⌶(S) ⊗ ⊥(S) == ⊥(S)#zero for multiplication
            end
        end
        for S in [TropicalCSemifield{T}, MinPlusCSemifield{T},ArcticCSemifield{T}, MaxPlusCSemifield{T}]
            for x in S.(-3.0:1.0:3.0)#finite elements
                #y = S(2.0)
                #z = S(3.0)
                #@test x ⊕ y == S(_logaddexp(τ, val(x), val(y)))
                #@test x ⊗ y == S(3.0)
                @test x ⊕ x == x#idempotency of addition
                #missing power.x ⊗ x == x^2
                @test x ⊗ ⌶(S) == x == ⌶(S) ⊗ x
                @test x ⊕ ⊤(S) == ⊤(S) == ⊤(S) ⊕ x
                @test x ⊗ ⊥(S) == ⊥(S) == ⊥(S) ⊗ x
            end
            @test ⊤(S) ⊕ ⊤(S) == ⊤(S)#idempotence of addition
            @test ⊤(S) ⊗ ⊤(S) == ⊤(S)#idempotence of multiplication
            @test ⊤(S) ⊕ ⊥(S) == ⊤(S)#null element for addition
            @test ⊤(S) ⊗ ⊥(S) == ⊥(S)#absorptive element for multiplication
            @test ⊥(S) ⊕ ⊤(S) == ⊤(S)# null element for addition, commutativity
            @test ⊥(S) ⊗ ⊤(S) == ⊥(S)#absorptive element for multiplication, commutativity
            @test ⊥(S) ⊕ ⌶(S) == ⌶(S)#identity element for addition
            @test ⊥(S) ⊗ ⌶(S) == ⊥(S)#identity element for multiplication
            @test ⌶(S) ⊕ ⊤(S) == ⊤(S)#absorptive element for addition
            @test ⌶(S) ⊗ ⊤(S) == ⊤(S)#identity element for multiplication
            @test ⌶(S) ⊕ ⊥(S) == ⌶(S)#identity element for addition
            @test ⌶(S) ⊗ ⊥(S) == ⊥(S)#zero for multiplication
        end
    end
end