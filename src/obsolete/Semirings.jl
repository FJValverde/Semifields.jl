# SPDX-License-Identifier: CECILL-2.1

module Semirings

using ChainRulesCore
import LogExpFunctions: logaddexp

export
    # Main API.
    ⊕,
    ⊗,
    val,

    # Partial derivative functions.
    ∂sum,
    ∂rmul,
    ∂lmul,

    # Concrete types.
    Semiring,
    ArcticSemiring,
    BoolSemiring,
    LogSemiring,
    ProbSemiring,
    TropicalSemiring


###############################################################################
# Abstract semiring type.

@doc raw"""
    abstract type Semiring end

Abstract type for a semiring ``(S, \oplus, \otimes, \bar{0}, \bar{1})``.
```math
(S, \\oplus, \\otimes, \\bar{0}, \\bar{1}
```

"""
abstract type Semiring{T} end

"""
    x ⊕ y

Semiring addition.
"""
⊕

"""
    x ⊗ y

Semiring multiplication.
"""
⊗

"""
    val(x::Semiring)

Return the "real" value / object wrapped in the semiring type
"""
val(x::Semiring) = x.val

"""
    Base.valtype(::Type{<:Semiring{T}}) where T

Return the type of the value wrapped by the semiring.
"""
Base.valtype(::Type{<:Semiring{T}}) where T = T

Base.zero(x::Semiring) = zero(typeof(x))
Base.one(x::Semiring) = one(typeof(x))

function Base.:*(i::Integer, s::Semiring)
    i < 0 && throw(ArgumentError("integer has to be positive or zero"))
    iszero(i) && return zero(s)

    res = zero(s)
    for n in 1:i
        res = res ⊕ s
    end
    res
end
Base.:*(s::Semiring, i::Integer) = i * s


Base.convert(T::Type{<:Semiring}, x::Number) = T(x)

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

⊕(x::BoolSemiring, y::BoolSemiring) = BoolSemiring(x.val || y.val)
⊗(x::BoolSemiring, y::BoolSemiring) = BoolSemiring(x.val && y.val)
Base.zero(::Type{<:BoolSemiring}) = BoolSemiring(false)
Base.one(::Type{<:BoolSemiring}) = BoolSemiring(true)

"""
    struct LogSemiring{T,τ} <: Semiring{T}
        val::T
    end

Logarithmic semiring: ``(\\mathbb{R} \\cup \\{- \\infty \\}, \\oplus_{\\log}, +, -\\infty, 0)``
where

```math
x \\oplus y = \\frac{1}{\\tau} \\log ( e^{\tau x} + e^{\tau y} ).
```
"""
struct LogSemiring{T<:AbstractFloat,τ} <: Semiring{T}
    val::T
end

_logaddexp(τ, x, y) = inv(τ) * logaddexp(τ*x, τ*y)

⊕(x::LogSemiring{T,τ}, y::LogSemiring{T,τ}) where {T,τ} = LogSemiring{T,τ}(_logaddexp(τ, val(x), val(y)))
⊗(x::S, y::S) where S<:LogSemiring = S(val(x) + val(y))
Base.zero(S::Type{<:LogSemiring{T,τ}}) where {T,τ} = S(ifelse(τ > 0, T(-Inf), T(Inf)))
Base.one(S::Type{<:LogSemiring}) = S(0)

∂sum(z::LogSemiring{T,τ}, x::LogSemiring{T,τ}) where {T,τ} =
    val(z) == -Inf ? zero(T) : exp(τ*(val(x) - val(z)))
∂rmul(x::S, a::S) where S<:LogSemiring = valtype(S)(1)
∂lmul(a::S, x::S) where S<:LogSemiring = valtype(S)(1)

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

∂sum(z::S, x::S) where S  = 1
∂rmul(x::S, a::S) where S<:ProbSemiring = val(a)
∂lmul(a::S, x::S) where S<:ProbSemiring = val(a)

"""
    const TropicalSemiring{T} = LogSemiring{T,-Inf} where T

Tropical semiring: ``(\\mathbb{R} \\cup \\{- \\infty \\}, min, +, \\infty, 0)``.
"""
const TropicalSemiring{T} = LogSemiring{T,-Inf} where T
⊕(x::S, y::S) where S<:TropicalSemiring = S(min(val(x), val(y)))

∂sum(z::S, x::S) where S<:TropicalSemiring = valtype(S)(x == z)

"""
    const ArcticSemiring{T} = LogSemiring{T,Inf} where T

Tropical semiring: ``\\langle \\mathbb{R} \\cup \\{\\infty \\}, max, +, -\\infty, 0 \\rangle``.
"""
const ArcticSemiring{T} = LogSemiring{T,Inf} where T
⊕(x::S, y::S) where S<:ArcticSemiring = S(max(val(x), val(y)))

∂sum(z::S, x::S) where S<:ArcticSemiring = valtype(S)(x == z)


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
        if x == y == -Inf
            (NoTangent(), NoTangent(), 0, 0)
        else
            (NoTangent(), NoTangent(), exp(τ*(x - z)) * z̄, exp(τ*(y - z)) * z̄)
        end
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
