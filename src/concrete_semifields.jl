import ..Semirings: ⊕,⊗,val, valtype,∂sum,∂rmul,∂lmul

###############################################################################
# Concrete semifield types.
"""
    struct BoolSemifield <: Semifield{Bool} <: Number
        val::Bool
    end

Boolean semifield: ``\\langle \\{0, 1\\}, \\lor, \\land, 0, 1 \\rangle``.
"""
struct BoolSemifield <: Semifield{Bool}
    val::Bool
end

⊕(x::BoolSemifield, y::BoolSemifield) = BoolSemifield(x.val || y.val)
⊗(x::BoolSemifield, y::BoolSemifield) = BoolSemifield(x.val && y.val)
Base.zero(::Type{<:BoolSemifield}) = BoolSemifield(false)
Base.one(::Type{<:BoolSemifield}) = BoolSemifield(true)
top(::Type{<:BoolSemifield}) = BoolSemifield(true)

Base.inv(x::Type{<:BoolSemifield}) = BoolSemifield(!x.val)

"""
    struct LogSemifield{T,τ} <: Semifield{T}
        val::T
    end

Logarithmic semifield: ``(\\mathbb{R} \\cup \\{- \\infty \\}, \\oplus_{\\log}, +, -\\infty, 0)``
where ``\\tau\\belongs\\overline\\mathbb{R}`` and ``\\oplus_{\\log}`` is the log-sum-exp operation:
```math
x \\oplus y = \\frac{1}{\\tau} \\log ( e^{\\tau x} + e^{\\tau y} ).
```
Note that ``\\tau \\neq 0'' and that ``\\tau < 0'' makes is a *cost* semifield, while
``\\tau > 0'' makes a *utility* semifield.
"""
struct LogSemifield{T<:AbstractFloat,τ} <: Semifield{T} 
    val::T
end

_logaddexp(τ, x, y) = inv(τ) * logaddexp(τ*x, τ*y)#operations in floats. 
#The inversion here is the reason why τ cannot be zero.

function ⊕(x::LogSemifield{T,τ}, y::LogSemifield{T,τ}) where {T,τ}
    iszero(x) && return y#early termination
    iszero(y) && return x#early termination
    #( x == zero(LogSemifield{T,τ}) ? y :
    #( y == zero(LogSemifield{T,τ}) ? x :
    return(LogSemifield{T,τ}(_logaddexp(τ, val(x), val(y))))
end


function ⊗(x::S, y::S) where S<:LogSemifield
    iszero(x) && return x#early termination
    iszero(y) && return y#early termination
    return S(val(x) + val(y))
    #(x == zero(S) ? x : (y == zero(S) ? y : S(val(x) + val(y))))
end

Base.inv(x::S) where S<:LogSemifield = S(-val(x))
#Base.inv(x::LogSemifield{T,τ}) where {T,τ}  = LogSemifield{T,τ}(-val(x))

function Base.zero(S::Type{<:LogSemifield{T,τ}}) where {T,τ}
    iszero(τ) && throw(ArgumentError("Rényi parameter has to be nonzero"))
    S(τ > 0.0 ? T(-Inf) : T(Inf))
end
#Base.zero(S::Type{<:LogSemifield{T,τ}}) where {T,τ} = S(ifelse(τ > 0, T(-Inf), T(Inf)))

Base.one(S::Type{<:LogSemifield{T,τ}}) where {T,τ}  = S(zero(T))

"""
    top(S::Type{<:LogSemifield{T,τ}}) where {T,τ}

Defined as the inverse of the zero.
"""
top(S::Type{<:LogSemifield{T,τ}}) where {T,τ} = S(τ > 0 ? T(Inf) : T(-Inf))

function ChainRulesCore.rrule(::typeof(_logaddexp), b, x, y)
    z = _logaddexp(b, x, y)
    
    function _logaddexp_pullback(z̄)#is this to avoid vanishing gradients?
    #FVA: the following function to avoid vanishing gradients implies \tau > 0
     #   if x == y == -Inf
     #       (NoTangent(), NoTangent(), 0, 0)
        if b > 0 && x == y == -Inf
            (NoTangent(), NoTangent(), 0, 0)
        elseif b < 0 && x == y == Inf
            (NoTangent(), NoTangent(), 0, 0)
        else
            (NoTangent(), NoTangent(), exp(b*(x - z)) * z̄, exp(b*(y - z)) * z̄)
        end
        #=
        if x == y == -Inf
            (NoTangent(), NoTangent(), 0, 0)
        else
            (NoTangent(), NoTangent(), exp(b*(x - z)) * z̄, exp(b*(y - z)) * z̄)
        end
        =#
    end
    z, _logaddexp_pullback
end
    
∂sum(z::LogSemifield{T,τ}, x::LogSemifield{T,τ}) where {T,τ} =
    #val(z) == -Inf ? zero(T) : exp(τ*(val(x) - val(z)))
    z == zero(LogSemifield{T,τ}) ? zero(T) : exp(τ*(val(x) - val(z)))
∂rmul(x::S, a::S) where S<:LogSemifield = valtype(S)(1)
∂lmul(a::S, x::S) where S<:LogSemifield = valtype(S)(1)

#=
# This probability semifield has problem. It is better to use the analog of max-min-times. 
"""
    struct ProbSemifield{T<:AbstractFloat} <: Semifield{T}
        val::T
    end

Probability semifield ``( (\\mathbb{R}_+``, +, \\cdot, 0, 1 )``.
"""
struct ProbSemifield{T<:AbstractFloat} <: Semifield{T}
    val::T
end

# ⊕(x::ProbSemifield, y::ProbSemifield) =
#     ( x == zero(ProbSemifield) ? y :
#     ( y == zero(ProbSemifield) ? x :
#     ProbSemifield(val(x) + val(y))))
⊕(x::S, y::S) where S<:ProbSemifield =
    ( x == zero(S) ? y :
    ( y == zero(S) ? x :
    S(val(x) + val(y))))

# ⊗(x::ProbSemifield, y::ProbSemifield) =
#     ( x == zero(ProbSemifield) ? x :
#     ( y == zero(ProbSemifield) ? y :
#     ProbSemifield(val(x) * val(y))))
⊗(x::S, y::S) where S<:ProbSemifield = S(val(x) * val(y))
    #iszero(x) || iszero(y) ? zero(S) : S(val(x) * val(y)) 
    #( x == zero(S) ? x :
    #( y == zero(S) ? y :
    #S(val(x) * val(y))))

# Watch out for this inversion, which might not be what you want, but the complementary, 
Base.inv(x::S) where S<:ProbSemifield = S(inv(val(x)))

Base.zero(S::Type{<:ProbSemifield{T}}) where T = S(zero(T))

Base.one(S::Type{<:ProbSemifield{T}}) where T = S(one(T))

top(S::Type{<:ProbSemifield}) = S(Inf)
#top(::ProbSemifield{T}) where T = inv(zero(T))#Two constructor calls

∂sum(z::S, x::S) where S  = 1
∂rmul(x::S, a::S) where S<:ProbSemifield = val(a)
∂lmul(a::S, x::S) where S<:ProbSemifield = val(a)
=#

"""
    const TropicalSemifield{T} = LogSemifield{T,-Inf} where T

Tropical semifield: ``(\\mathbb{R} \\cup \\{- \\infty \\}, min, +, \\infty, 0)``.

This is the Logexp semifielf with ``\\tau = -\\infty``.
"""
const TropicalSemifield{T} = LogSemifield{T,-Inf} where T

⊕(x::S, y::S) where S<:TropicalSemifield = S(min(val(x), val(y)))

∂sum(z::S, x::S) where S<:TropicalSemifield = valtype(S)(x == z)

"""
    const ArcticSemifield{T} = LogSemifield{T,Inf} where T

Tropical semifield: ``\\langle \\mathbb{R} \\cup \\{\\infty \\}, max, +, -\\infty, 0 \\rangle``.
This is the Logexp semifielf with ``\\tau = \\infty``.

"""
const ArcticSemifield{T} = LogSemifield{T,Inf} where T
⊕(x::S, y::S) where S<:ArcticSemifield = S(max(val(x), val(y)))

# WRONG after the original author of Semirings. 
∂sum(z::S, x::S) where S<:ArcticSemifield = valtype(S)(x == z)

"""
    MaxMinTimesSemifield{T} where T
"""
const MaxTimesSemifield{T} = LogSemifield{T,Inf} where T
"""
    MinTimesSemifield{T} where T
"""
const MinTimesSemifield{T} = LogSemifield{T,-Inf} where T