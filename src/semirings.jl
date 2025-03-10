# SPDX-License-Identifier: CECILL-2.1
using Reexport
@reexport module Semirings
"""
    Semirings   

A module to define generic and concrete Semirings.
"""

#@reexport module <modulename> ... end defines module <modulename> 
# and also re-exports its symbols:
using ChainRulesCore
import LogExpFunctions: logaddexp

export
    # Main API for semirings.
    ⊕,
    ⊗,
    val, valtype,

    # Partial derivative functions.
    ∂sum,
    ∂rmul,
    ∂lmul,

    # Concrete types.
    Semiring,
    BoolSemiring,
    EntropySemiring,
    #ProbSemiring,
    EnergySemiring,
    ArcticSemiring,
    TropicalSemiring


###############################################################################
# Abstract semiring type.

"""
    abstract type Semiring end

Abstract type for a semiring ``(S, \\oplus, \\otimes, \\bar{0}, \\bar{1})``.
```math
(S, \\oplus, \\otimes, \\bar{0}, \\bar{1})
```
"""
abstract type Semiring{T} end

"""
    x ⊕ y

An abstract semiring addition. ``\\langle S, ⊕\\rangle`` 
is typically a commutative monoid. 
"""
⊕

"""
    x ⊗ y

An abstract semiring addition. ``\\langle S, ⊗\\rangle`` 
is typically a commutative monoid. 
"""
⊗

"""
    val(x::Semiring{T}) where T → T 

Return the "real" value / object wrapped in the semiring type. 
"""
val(x::Semiring) = x.val

"""
    Base.valtype(::Type{<:Semiring{T}}) where T → T 

Return the type of the value wrapped by the semiring.
"""
Base.valtype(::Type{<:Semiring{T}}) where T = T

"""
    Base.zero(::Type{<:Semiring{T}}) where T → Semiring{T}

The additive zero element in the semiring. 
"""
Base.zero(x::Semiring) = zero(typeof(x))
"""
    Base.one(::Type{<:Semiring{T}}) where T → Semiring{T}

The multiplicative unit element in the semiring. 
""" 
Base.one(x::Semiring) = one(typeof(x))

"""
    Base.:⊗(i::Integer, s::Semiring) → Semiring

A commutative multiplication by integers, not really know where it leads to.
"""
#function Base.:⊗(i::Integer, s::Semiring)
#function Base.:⊗(i::T, s::Semiring{T}) where T
    #i < 0 && throw(ArgumentError("integer has to be positive or zero"))
    #iszero(i) && return zero(s)

    #res = zero(s)
    #for n in 1:i
    #    res = res ⊕ s
    #end
    #res
#end
#Base.:⊗(s::Semiring, i::Integer) = i ⊗ s


"""
    Base.convert(T::Type{<:Semiring}, x::Number)

TBW
"""
Base.convert(T::Type{<:Semiring}, x::Number) = T(x)
Base.convert(T::Type{<:Semiring}, x::Semiring) = T(x.val)

Base.show(io::IO, x::Semiring) = print(io, val(x))


###############################################################################
# Concrete semiring types.

"""
    struct BoolSemiring <: Semiring{Bool}
        val::Bool
    end

Boolean semiring: ``\\langle \\{0, 1\\}, \\lor, \\land, 0, 1 \\rangle``.
"""
struct BoolSemiring <: Semiring{Bool}
    val::Bool
end

"""
    x ⊕ y -> BoolSemiring

Or operation in the Boolean semiring.
"""
⊕(x::BoolSemiring, y::BoolSemiring) = BoolSemiring(x.val || y.val)

"""
    x ⊗ y -> BoolSemiring

And operation in The Boolean semiring.
"""
⊗(x::BoolSemiring, y::BoolSemiring) = BoolSemiring(x.val && y.val)

"""
    Base.zero(::Type{<:BoolSemiring}) → BoolSemiring

The zero bit is also the additive zero.
"""
Base.zero(::Type{<:BoolSemiring}) = BoolSemiring(false)

"""
    Base.one(::Type{<:BoolSemiring}) → BoolSemiring

The one bit is also the multiplicative unit.
"""
Base.one(::Type{<:BoolSemiring}) = BoolSemiring(true)

"""
    struct EntropySemiring{T,τ} <: Semiring{T}
        val::T
    end

Entropic semiring: ``(\\mathbb{R} \\cup \\{- \\infty \\}, \\oplus_{\\tau}, +, -\\infty, 0)``
where
```math
x \\oplus_{\\tau} y = \\frac{1}{\\tau} \\log ( e^{\\tau x} + e^{\\tau y} ).
```
"""
struct EntropySemiring{T<:AbstractFloat,τ} <: Semiring{T}
    val::T
end

_logaddexp(τ, x, y) = inv(τ) * logaddexp(τ*x, τ*y)

"""
    x ⊕ y -> EntropySemiring{T,τ}

Information addition in the entropy semiring. τ is the Rényi parameter.
"""
⊕(x::EntropySemiring{T,τ}, y::EntropySemiring{T,τ}) where {T,τ} = EntropySemiring{T,τ}(_logaddexp(τ, val(x), val(y)))
"""
    x ⊗ y -> EntropySemiring{T,τ}

Information multiplication in the entropy semiring. Does not depend on τ.
"""
⊗(x::EntropySemiring{T,τ}, y::EntropySemiring{T,τ}) where {T,τ} = EntropySemiring{T,τ}(val(x) + val(y))
⊗(x::S, y::S) where S<:EntropySemiring = S(val(x) + val(y))

"""
    Base.zero(::Type{<:EntropySemiring{T,τ}}) where {T,τ} → EntropySemiring{T,τ}
"""
Base.zero(S::Type{<:EntropySemiring{T,τ}}) where {T,τ} = S(ifelse(τ > 0, T(-Inf), T(Inf)))

"""
    Base.one(::Type{<:EntropySemiring{T}}) where T → EntropySemiring{T}

The multiplicative unit in the entropy semiring.
"""
Base.one(S::Type{<:EntropySemiring}) = S(0)

∂sum(z::EntropySemiring{T,τ}, x::EntropySemiring{T,τ}) where {T,τ} =
    val(z) == -Inf ? zero(T) : exp(τ*(val(x) - val(z)))
∂rmul(x::S, a::S) where S<:EntropySemiring = valtype(S)(1)
∂lmul(a::S, x::S) where S<:EntropySemiring = valtype(S)(1)

"""
    const TropicalSemiring{T} = EntropySemiring{T,-Inf} where T

Tropical semiring: ``(\\mathbb{R} \\cup \\{- \\infty \\}, min, +, \\infty, 0)``.
"""
const TropicalSemiring{T} = EntropySemiring{T,-Inf} where T
⊕(x::S, y::S) where S<:TropicalSemiring = S(min(val(x), val(y)))

#=
∂sum(z::S, x::S) where S<:TropicalSemiring = valtype(S)(x == z)
=#
"""
    const ArcticSemiring{T} = EntropySemiring{T,Inf} where T

Tropical semiring: ``\\langle \\mathbb{R} \\cup \\{\\infty \\}, max, +, -\\infty, 0 \\rangle``.
"""
const ArcticSemiring{T} = EntropySemiring{T,Inf} where T
⊕(x::S, y::S) where S<:ArcticSemiring = S(max(val(x), val(y)))

#=
∂sum(z::S, x::S) where S<:ArcticSemiring = valtype(S)(x == z)
=#

"""
    struct EnergySemiring{T<:AbstractFloat} <: Semiring{T}
        val::T
    end

Energy semiring ``( (\\mathbb{R}_+``, +, \\cdot, 0, 1 )``.
"""
struct EnergySemiring{T<:AbstractFloat} <: Semiring{T}
    val::T
end

⊕(x::EnergySemiring, y::EnergySemiring) = EnergySemiring(val(x) + val(y))
⊗(x::EnergySemiring, y::EnergySemiring) = EnergySemiring(val(x) * val(y))
Base.zero(S::Type{<:EnergySemiring{T}}) where T = S(zero(T))
Base.one(S::Type{<:EnergySemiring{T}}) where T = S(one(T))

∂sum(z::S, x::S) where S  = 1#one(S)
∂rmul(x::S, a::S) where S<:EnergySemiring = val(a)
∂lmul(a::S, x::S) where S<:EnergySemiring = val(a)

#=
#Energy semirings are missing here.

"""
    struct ProbSemiring{T<:AbstractFloat} <: Semiring{T}
        val::T
    end

Probability semiring ``( (\\mathbb{R}_+``, +, \\cdot, 0, 1 )``.
"""
struct ProbSemiring{T<:AbstractFloat} <: Semiring{T}
    val::T
end

⊕(x::ProbSemiring, y::ProbSemiring) = ProbSemiring(val(x) + val(y))
⊗(x::ProbSemiring, y::ProbSemiring) = ProbSemiring(val(x) * val(y))
Base.zero(S::Type{<:ProbSemiring{T}}) where T = S(zero(T))
Base.one(S::Type{<:ProbSemiring{T}}) where T = S(one(T))
=#
#=
∂sum(z::S, x::S) where S  = 1
∂rmul(x::S, a::S) where S<:ProbSemiring = val(a)
∂lmul(a::S, x::S) where S<:ProbSemiring = val(a)
=#

###############################################################################
# Differentiation.

"""
    ∂sum(z, x)

Compute the partial derivative of `z = x ⊕ y ⊕ z ⊕ ...` w.r.t. `x`.
"""
∂sum

"""
    ∂rmul(x, a)

Compute the partial derivative of `x ⊗ a` w.r.t. `x`.
"""
∂rmul

"""
    ∂lmul(a, x)

Compute the partial derivative of `a ⊗ x` w.r.t. `x`.
"""
∂lmul

function ChainRulesCore.rrule(::typeof(_logaddexp), τ, x, y)
    z = _logaddexp(τ, x, y)

    function _logaddexp_pullback(z̄)
        #FVA: the following has to take into consideration the value
        # of the Rényi parameter τ.
        if τ > 0 && x == y == -Inf
            (NoTangent(), NoTangent(), 0, 0)
        elseif τ < 0 && x == y == Inf
            (NoTangent(), NoTangent(), 0, 0)
        else
            (NoTangent(), NoTangent(), exp(b*(x - z)) * z̄, exp(b*(y - z)) * z̄)
        end
        #=
        if x == y == -Inf
            (NoTangent(), NoTangent(), 0, 0)
        else
            (NoTangent(), NoTangent(), exp(τ*(x - z)) * z̄, exp(τ*(y - z)) * z̄)
        end
        =#
    end

    z, _logaddexp_pullback
end

function ChainRulesCore.rrule(::typeof(val), a::Semiring)
    b = val(a)
    function val_pullback(b̄)
        ZeroTangent(), Tangent{typeof(a)}(; val=b̄)
    end
    b, val_pullback
end

function ChainRulesCore.rrule(S::Type{<:Semiring}, a)
    b = S(a)
    function semiring_pullback(b̄)
        ZeroTangent(), b̄.val
    end
    b, semiring_pullback
end

end
