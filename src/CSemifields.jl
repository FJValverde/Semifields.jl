#FVA: once CSemifields are working, since the focus of the project is on double complete
# semifields, the following module should be re-exported.
#
#using Reexport
#@reexport 
module CSemifields
"""
    CSemifields -- Complete Semifields

A module to define complete semifields as pairs of two semifields
related by inversion.
"""

export top, # the top element of the complete semifield.
    ⊤,# an alternate name for the top element.
    ⊥, #an alternate name for the bottom element. 
    ⌶,#\topbot, an alternative for the unit element

    CSemifield,#the basic abstract type of Complete(d) semifields.
    TernaryCSemifield,#The basic semifield included in all the rest.
    EntropyCSemifield,#The pair of semifields of logarithmic ranges.
    TropicalCSemifield, MinPlusCSemifield,#The pair of semifields of min-plus ranges.
    ArcticCSemifield, MaxPlusCSemifield#The pair of semifields of max-plus ranges.
#=
    MaxMinTimes
=#
# CAVEAT! The complete Semifields are never subtypes of the 
# concrete semifields

#using ChainRulesCore
#import LogExpFunctions: logaddexp
#using Reexport
#@reexport using ..Semifields#This should not be needed, since Semirings is re-exported by Semifields.

import ..Semirings: ⊕,⊗,val, valtype,∂sum,∂rmul,∂lmul, Semiring, _logaddexp, _holderaverage

#import ..GenericSemifields: Semifield
import Semifields: Semifield, EntropySemifield

"""
    CSemifield{T} <: Semifield{T}

Complete semifields with a top added. 
"""
abstract type CSemifield{T} <: Semifield{T} 
end

"""
    ⊤(x::CSemifield) → CSemifield
    ⊥(x::CSemifield) → CSemifield
    ⌶(x::CSemifield) → CSemifield  

Visual notation for top(⊤, \\top), bottom(⊥, \\bot) and unit (⌶, \\topbot) elements of the semifield.
"""
⊤(S::Type{<:CSemifield}) = top(S)
⊥(S::Type{<:CSemifield}) = zero(S)
⌶(S::Type{<:CSemifield}) = one(S)

"""
    convert(S::Type{<:CSemifield}, x::Number) → CSemifield

Basic convert function for Semifields elements.
"""
Base.convert(S::Type{<:CSemifield}, x::Number) = S(x)#For Number types
Base.convert(S::Type{<:CSemifield}, x::CSemifield) = S(val(x))#For subtypes of semifields.


#=
"""
    TernarySemifield <: CSemifield{Some(Bool)}

This is not quite the same as the BoolSemifield, that is already a complete 
semifield. It is enriched with a "middle" element that is not manifested, actually, 
as showing "no polarity" in the ternary logic. 
"""
struct TernarySemifield <: CSemifield{Some{Bool}}
    val::Some{Bool}
end
#Base.zero(::Some{Bool}) = TernarySemifield(false)
Base.zero(::Type{<:TernarySemifield}) = TernarySemifield(Some(false))
top(::Type{<:TernarySemifield}) = TernarySemifield(Some(true))
Base.one(::Type{<:TernarySemifield}) = TernarySemifield(nothing)
Base.inv(x::TernarySemifield) = isnothing(x.val) ? x : TernarySemifield(!x.val)

=#
"""
    TernaryCSemifield <: CSemifield{Union{Bool,Nothing}}

This is not quite the same as the BoolSemifield, that is already a complete 
semifield. It is enriched with a "middle" element that is not manifested, actually, 
as showing "no polarity" in the ternary logic. It is ordered from bottom to unit to top. 
"""
struct TernaryCSemifield <: CSemifield{Union{Bool,Nothing}}
    val::Union{Bool,Nothing}
end

Base.zero(S::Type{<:TernaryCSemifield}) = S(false)
top(S::Type{<:TernaryCSemifield}) = S(true)
Base.one(S::Type{<:TernaryCSemifield}) = S(nothing)
Base.inv(x::TernaryCSemifield) = isnothing(x.val) ? x : TernaryCSemifield(!x.val)

#=
"""
    ⊤(x::TernaryCSemifield) → TernaryCSemifield

An alias for the top element of the TernaryCSemifield.
"""
⊤(x::TernaryCSemifield) = top(typeof(x))

"""
    ⊥(x::TernaryCSemifield) → TernaryCSemifield

An alias for the zero element of the TernaryCSemifield.
"""
⊥(x::TernaryCSemifield) = zero(typeof(x))
=#
⊕(x::TernaryCSemifield, y::TernaryCSemifield) = 
    x.val == false ? y : 
    y.val == false ? x :
    x.val == true ? x :
    y.val == true ? y : x# == y == nothing

⊗(x::TernaryCSemifield, y::TernaryCSemifield) = 
    x.val == false ? x : 
    y.val == false ? y :
    isnothing(x.val) ? y :
    isnothing(y.val) ? x : nothing

#FVA: the doc raw prevents the docstrings from being interpreted.  

@doc raw"""
struct EntropyCSemifield{T,τ} <: CSemifield{T}
    val::T
end

Entropic complete semifields: ``(\\mathbb{R} \\cup \\{- \\infty \\}, \\oplus_{\\log}, +, -\\infty, 0, \\infty)``
where ``\tau\belongsto\overline\mathbb{R}`` and ``\oplus_{\log}`` is the log-sum-exp operation:
```math
x \oplus y = \frac{1}{\tau} \log ( e^{\tau x} + e^{\tau y} ).
```
Note that ``\tau \neq 0`` and that ``\tau < 0`` makes is a *cost* semifield, while
``\tau > 0`` makes a *utility* semifield.
"""
struct EntropyCSemifield{T,τ} <: CSemifield{T} 
    val::T
end

"""
    zero(S::Type{<:EntropyCSemifield{T,τ}}) → EntropyCSemifield{T,τ}

The zero in the semifield, when τ != 0.0.
"""
Base.zero(S::Type{<:EntropyCSemifield{T,τ}}) where {T,τ} = S(-sign(τ)*typemax(T))
#=
function Base.zero(S::Type{<:EntropyCSemifield{T,τ}}) where {T,τ}
    S( 
        τ > 0.0 ? T(-Inf) : 
        τ < 0.0 ? T(Inf) : 
        throw(ArgumentError("Rényi parameter has to be nonzero"))
    )
end
=#
#Base.zero(S::Type{<:EntropySemifield{T,τ}}) where {T,τ} = S(ifelse(τ > 0, T(-Inf), T(Inf)))

"""
    one(S::Type{<:EntropyCSemifield{T,τ}}) → EntropyCSemifield{T,τ}

The unit in the semifield, when τ != 0.0.
"""
Base.one(S::Type{<:EntropyCSemifield{T,τ}}) where {T,τ} = S(zero(T))

"""
    top(S::Type{<:EntropyCSemifield{T,τ}}) where {T,τ}

Defined as the inverse of the zero.
```julia
using Tests
K = EntropyCSemifield{Float64,1.0}
@test ⊤(K) == inv(zero(K))#true
@test ⊤(K) == top(K)
@test ⊤(K) == inv(⊥(K))
```
"""
top(S::Type{<:EntropyCSemifield{T,τ}}) where {T,τ} = S(τ > 0 ? T(Inf) : T(-Inf))

"""
    inv(x::EntropyCSemifield{T,τ}) where {T,τ}  = EntropyCSemifield{T,τ}(-val(x))

The multiplicative inverse unary operator. 
"""
Base.inv(S::Type{<:EntropyCSemifield{T,τ}}) where {T,τ}  = S(-val(x))

"""
    ⊕(x::EntropyCSemifield{T,τ}, y::EntropyCSemifield{T,τ}) where {T,τ} = 

Addition in complete semifields, including results for the top which is absorptive. 
"""
function ⊕(x::EntropyCSemifield{T,τ}, y::EntropyCSemifield{T,τ}) where {T,τ}
    #iszero(x) && return y#early termination
    #iszero(y) && return x#early termination
    #return(EntropyCSemifield{T,τ}(_logaddexp(τ, val(x), val(y))))
    x == zero(x) ? y :
    y == zero(y) ? x :
    EntropyCSemifield{T,τ}(_logaddexp(τ, val(x), val(y)))
end

"""
    ⊗(x::S, y::S) where S<:EntropySemifield
"""
function ⊗(x::S, y::S) where S<:EntropyCSemifield
    #iszero(x) && return x#early termination
    #iszero(y) && return y#early termination
    #return S(val(x) + val(y))
    x == zero(S) ? x : 
    y == zero(S) ? y : 
    S(val(x) + val(y))
end #module 

"""
    const TropicalCSemifield{T} = EntropyCSemifield{T,-Inf} where T

Completed tropical or min-plus semifield: ``(\\mathbb{R} \\cup \\{- \\infty \\}, min, +, \\infty, 0, -\\infty)``.

This is the completed entropy semifield with ``\\tau = -\\infty``.
"""
const TropicalCSemifield{T} = EntropyCSemifield{T,-Inf} where T
const MinPlusCSemifield{T} = TropicalCSemifield{T} where T

⊕(x::S, y::S) where S<:TropicalCSemifield = S(min(val(x), val(y)))
#⊗(x::S, y::S) where S<:TropicalCSemifield = S(val(x) + val(y))#Not changed from the original
Base.zero(S::Type{<:TropicalCSemifield{T}}) where T = S(Inf)#Shorter than the original
#Base.one(S::Type{<:TropicalCSemifield{T}}) where T = S(zero(T))#Not changed from the original
top(S::Type{<:TropicalCSemifield{T}}) where T = S(-Inf)#Shorter than the original
#Base.inv(x::S) where S<:TropicalCSemifield = S(-val(x))
#∂sum(z::S, x::S) where S<:TropicalCSemifield = valtype(S)(x == z)
#∂rmul(x::S, a::S) where S<:TropicalCSemifield = valtype(S)(1)#FVA: not sure of this
#∂lmul(a::S, x::S) where S<:TropicalCSemifield = valtype(S)(1)#FVA: not sure of this

"""
    const ArcticCSemifield{T} = EntropyCSemifield{T,Inf} where T

Complete arctic or max-plus semifield: ``\\langle \\mathbb{R} \\cup \\{\\infty \\}, max, +, -\\infty, 0, +\\infty \\rangle``.
This is the complete entropy semifield with ``\\tau = \\infty``.

"""
const ArcticCSemifield{T} = EntropyCSemifield{T,Inf} where T#Technically correct, but a mouthful
const MaxPlusCSemifield{T} = ArcticCSemifield{T} where T#Preferred named
⊕(x::S, y::S) where S<:ArcticCSemifield = S(max(val(x), val(y)))
#⊗(x::S, y::S) where S<:ArcticCSemifield = S(val(x) + val(y))#Not changed from the original
Base.zero(S::Type{<:ArcticCSemifield{T}}) where T = S(-Inf)#Shorter than the original
#Base.one(S::Type{<:ArcticCSemifield{T}}) where T = S(zero(T))#Not changed from the original
top(S::Type{<:ArcticCSemifield{T}}) where T = S(Inf)#Shorter than the original
#Base.inv(x::S) where S<:ArcticCSemifield = S(-val(x))#Not changed from the original
#∂sum(z::S, x::S) where S<:ArcticCSemifield = valtype(S)(x == z)#FVA:  WRONG after the original author of Semirings. 
#∂rmul(x::S, a::S) where S<:ArcticCSemifield = valtype(S)(1)#FVA: not sure of this
#∂lmul(a::S, x::S) where S<:ArcticCSemifield = valtype(S)(1)#FVA: not sure of this


end#CSemifields