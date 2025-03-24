using Reexport
@reexport module DoubleCSemifields
"""
    DoubleCSemifields -- DoubleComplete Semifields

``(S, \\oplus, \\otimes, inv, zero, one, \\top)``

A module to define complete semifields as pairs of two semifields related by inversion. 

Double semifields, as use, e.g in convex analysis, max-min-plus algebra, probability theory, etc.
"""

export DoubleCSemifield,#the basic abstract type of Complete(d) semifields.
    #top, # the top element of the complete semifield.
    #⊤,# an alternate name for the top element.
    #⊥, #an alternate name for the bottom element. 
    #⌶,#\topbot, an alternative for the unit element
    #operations
    ⊕̇, ⊗̇,#upper addition ''\oplus\dot'', multiplication \otimes\dot 
    ⊕̧, ⊗̧,#lower addition ''\oplus\c'', multplication \otimes\c (for the time being)
    # TODO: try to actually find something like \oplus\ldot lower dot, instead of \c

    #TernaryDoubleCSemifield,#The basic semifield included in all the rest.
    EntropyDoubleCSemifield,#The pair of semifields of logarithmic ranges.
    MaxMinPlusSemifield#The pair of semifields of min-plus ranges.
    #MinMaxPlusSemifield#The pair of semifields of max-plus ranges.
#=
    EnergyDoubleCSemifield,#The pair of semifields of energy and power ranges.
    MaxMinTimesCSemifield
=#
#using ChainRulesCore
#import LogExpFunctions: logaddexp
using Reexport
#@reexport using ..Semifields#This should not be needed, since Semirings is re-exported by Semifields.

@reexport import ..Semirings: ⊕,⊗,val, valtype,∂sum,∂rmul,∂lmul, Semiring, _logaddexp, _holderaverage

@reexport import Semifields: Semifield, EntropySemifield, CSemifields, inv, zero, one

@reexport import ..CSemifields: CSemifield, EntropyCSemifield, top, ⊤, ⊥, ⌶

"""
    DoubleCSemifield{T} <: CSemifield{T}

Double complete semifields are amalgamated pairs of two semifields related by inversion, with
signature: ``(S, \\oplus\\c, \\oplus\\dot, \\otimes\\c, \\otimes\\dot, inv, zero, one, \\top)``

We actually need two complete semifields on the same carrier set but opposite order to build 
them. They have:
- A top element, ⊤, that is the inverse of the bottom element of the other semifield.
- A bottom element, ⊥, that is the inverse of the top element of the other semifield.
- A unit element, ⌶, that is the same in both semifields.
- An inversion operator, inv, that is the same in both semifields.
- ⊕̧ ``\\oplus\\c`` a lower addition operator, whose identity is the bottom element of the semifield.
- ⊗̧ ``\\otimes\\c`` a lower multiplication operator, for which the bottom is the zero element.
- ⊕̇ ``\\oplus\\dot`` an upper addition operator, whose identity is the bottom element.
- ⊗̇ ``\\otimes\\dot`` an upper multiplication operator, for which the top is the zero element.
"""
abstract type DoubleCSemifield{T} <: CSemifield{T} 
end

"""
    ⊤(x::DoubleCSemifield) → DoubleCSemifield
    ⊥(x::DoubleCSemifield) → DoubleCSemifield
    ⌶(x::DoubleCSemifield) → DoubleCSemifield  

Visual notation for top(⊤, \\top), bottom(⊥, \\bot) and unit (⌶, \\topbot) elements of the semifield.
"""
⊤(S::Type{<:DoubleCSemifield}) = top(S)
⊥(S::Type{<:DoubleCSemifield}) = zero(S)
⌶(S::Type{<:DoubleCSemifield}) = one(S)

"""
    convert(S::Type{<:DoubleCSemifield}, x::Number) → DoubleCSemifield

Basic convert function for Semifields elements.
"""
Base.convert(S::Type{<:DoubleCSemifield}, x::Number) = S(x)#For Number types
Base.convert(S::Type{<:DoubleCSemifield}, x::DoubleCSemifield) = S(val(x))#For subtypes of semifields.


#=
"""
    TernarySemifield <: DoubleCSemifield{Some(Bool)}

This is not quite the same as the BoolSemifield, that is already a complete 
semifield. It is enriched with a "middle" element that is not manifested, actually, 
as showing "no polarity" in the ternary logic. 
"""
struct TernarySemifield <: DoubleCSemifield{Some{Bool}}
    val::Some{Bool}
end
#Base.zero(::Some{Bool}) = TernarySemifield(false)
Base.zero(::Type{<:TernarySemifield}) = TernarySemifield(Some(false))
top(::Type{<:TernarySemifield}) = TernarySemifield(Some(true))
Base.one(::Type{<:TernarySemifield}) = TernarySemifield(nothing)
Base.inv(x::TernarySemifield) = isnothing(x.val) ? x : TernarySemifield(!x.val)


"""
    TernaryDoubleCSemifield <: DoubleCSemifield{Union{Bool,Nothing}}

This is not quite the same as the BoolSemifield, that is already a complete 
semifield. It is enriched with a "middle" element that is not manifested, actually, 
as showing "no polarity" in the ternary logic. It is ordered from bottom to unit to top.

There is no dual of this guy, as yet. It is order-aligned with the usual BoolSemifield.
"""
struct TernaryDoubleCSemifield <: DoubleCSemifield{Union{Bool,Nothing}}
    val::Union{Bool,Nothing}
end

Base.zero(S::Type{<:TernaryDoubleCSemifield}) = S(false)
top(S::Type{<:TernaryDoubleCSemifield}) = S(true)
Base.one(S::Type{<:TernaryDoubleCSemifield}) = S(nothing)
Base.inv(x::TernaryDoubleCSemifield) = isnothing(x.val) ? x : TernaryDoubleCSemifield(!x.val)

"""
    ⊕̇(x::TernaryDoubleCSemifield, y::TernaryDoubleCSemifield)

Lower addition, aka 'or' operation.
"""
⊕̧(x::TernaryDoubleCSemifield, y::TernaryDoubleCSemifield) = 
    x.val == false ? y : 
    y.val == false ? x :
    x.val == true ? x :
    y.val == true ? y : x# == y == nothing

⊗(x::TernaryDoubleCSemifield, y::TernaryDoubleCSemifield) = 
    x.val == false ? x : 
    y.val == false ? y :
    isnothing(x.val) ? y :
    isnothing(y.val) ? x : nothing
=#

"""
    struct EntropyDoubleCSemifield{T,τ} <: DoubleCSemifield{T}
        val::T
    end

Entropic double complete semifields: ``(\\mathbb{R} \\cup \\{- \\infty \\}, \\oplus_{\\log}, \\otimes_{\\log}, \\oplus^{\\log}, \\otimes^{\\log}, -\\infty, 0, \\infty)``
where ``\\tau\\belongs\\overline\\mathbb{R}`` and ``\\oplus_{\\log}`` is the log-sum-exp operation:
```math
x \\oplus_{\\tau} y = \\frac{1}{\\tau} \\log ( e^{\\tau x} + e^{\\tau y} ).
```
Note that ``\\tau \\neq 0'' and that ``\\tau < 0'' makes is a *cost* semifield, while
``\\tau > 0'' makes a *utility* semifield.
"""
struct EntropyDoubleCSemifield{T,τ} <: DoubleCSemifield{T} 
    val::T
end

#Base.convert(S::Type{<:EntropyDoubleCSemifield}, x::Number) = S(x)
#Base.convert(S::Type{<:EntropyDoubleCSemifield}, x::EntropyCSemifield) = S(val(x))#For subtypes of semifields.

"""
    zero(S::Type{<:EntropyDoubleCSemifield{T,τ}}) → EntropyDoubleCSemifield{T,τ}

The zero in the semifield, when τ != 0.0.
"""
Base.zero(S::Type{<:EntropyDoubleCSemifield{T,τ}}) where {T,τ} = S(-sign(τ)*typemax(T))
#=
function Base.zero(S::Type{<:EntropyDoubleCSemifield{T,τ}}) where {T,τ}
    S( 
        τ > 0.0 ? T(-Inf) : 
        τ < 0.0 ? T(Inf) : 
        throw(ArgumentError("Rényi parameter has to be nonzero"))
    )
end
=#
#Base.zero(S::Type{<:EntropySemifield{T,τ}}) where {T,τ} = S(ifelse(τ > 0, T(-Inf), T(Inf)))

"""
    one(S::Type{<:EntropyDoubleCSemifield{T,τ}}) → EntropyDoubleCSemifield{T,τ}

The multiplicative unit in the double semifield.
"""
Base.one(S::Type{<:EntropyDoubleCSemifield{T,τ}}) where {T,τ}  = S(zero(T))

"""
    top(S::Type{<:EntropyDoubleCSemifield{T,τ}}) where {T,τ}

Defined as the inverse of the zero.
```julia
using Tests
K = EntropyDoubleCSemifield{Float64,1.0}
@test ⊤(K) == inv(zero(K))#true
@test ⊤(K) == top(K)
@test ⊤(K) == inv(⊥(K))
```
"""
CSemifields.top(S::Type{<:EntropyDoubleCSemifield{T,τ}}) where {T,τ} = S(τ > 0 ? T(Inf) : T(-Inf))

"""
    inv(x::EntropyDoubleCSemifield{T,τ}) where {T,τ}  = EntropyDoubleCSemifield{T,τ}(-val(x))

The multiplicative inverse unary operator.
"""
Base.inv(x::S) where S <: EntropyDoubleCSemifield  = S(-val(x))


"""
    ⊗(x::EntropyDoubleCSemifield{T,τ}, y::EntropyDoubleCSemifield{T,τ}) where {T,τ}

Normal addition in the completed double semifield is the same ad the lower addition.
"""
function ⊕(x::EntropyDoubleCSemifield{T,τ}, y::EntropyDoubleCSemifield{T,τ}) where {T,τ}
    #iszero(x) && return y#early termination
    #iszero(y) && return x#early termination
    #return(EntropyDoubleCSemifield{T,τ}(_logaddexp(τ, val(x), val(y))))
    x == zero(x) ? y :
    y == zero(y) ? x :
    EntropyDoubleCSemifield{T,τ}(_logaddexp(τ, val(x), val(y)))
end

"""
    ⊕̧(x::EntropyDoubleCSemifield{T,τ}, y::EntropyDoubleCSemifield{T,τ}) where {T,τ}

Lower addition in double complete semifields, including results for the top which is absorptive. 
"""
⊕̧(x::EntropyDoubleCSemifield{T,τ}, y::EntropyDoubleCSemifield{T,τ}) where {T,τ} = x ⊕ y


"""
    ⊗(x::S, y::S) where S<:EntropySemifield

Lower multiplication in double complete semifields, including results for the top which is absorptive.
"""
function ⊗(x::S, y::S) where S<:EntropyDoubleCSemifield
    #iszero(x) && return x#early termination
    #iszero(y) && return y#early termination
    #return S(val(x) + val(y))
    x == zero(S) ? x : 
    y == zero(S) ? y : 
    S(val(x) + val(y))
end #module 

"""
    ⊗̧(x::EntropyDoubleCSemifield{T,τ}, y::EntropyDoubleCSemifield{T,τ}) where {T,τ}

Lower multiplication in double complete semifields, including results for the top which is absorptive.
"""
⊗̧(x::EntropyDoubleCSemifield{T,τ}, y::EntropyDoubleCSemifield{T,τ}) where {T,τ} = x ⊗ y


"""
    ⊕̇(x::EntropyDoubleCSemifield{T,τ}, y::EntropyDoubleCSemifield{T,τ}) where {T,τ}

Upper addition in double complete semifields, including results for the top which is the zero for it.
"""
function ⊕̇(x::S, y::S) where S<:EntropyDoubleCSemifield{T,τ} where {T,τ}
    #iszero(x) && return y#early termination
    #iszero(y) && return x#early termination
    #return(EntropyDoubleCSemifield{T,τ}(_logaddexp(τ, val(x), val(y))))
    x == top(S) ? y :
    y == top(S) ? x :
    #EntropyDoubleCSemifield{T,τ}(_logaddexp(-τ, val(x), val(y)))
    S(_logaddexp(-τ, val(x), val(y)))#This is the crux of the doubling!
end

"""
    oplusu(x::EntropyDoubleCSemifield{T,τ}, y::EntropyDoubleCSemifield{T,τ}) where {T,τ}

The upper multiplication.
"""
function oplusu(x::S, y::S) where S<:EntropyDoubleCSemifield
#function ⊗̇(x::S, y::S) where S <: EntropyCSemifield{T,τ} where {T,τ}
    #iszero(x) && return x#early termination
    #iszero(y) && return y#early termination
    #return S(val(x) + val(y))
    x == top(S) ? x : 
    y == top(S) ? y : 
    S(val(x) + val(y))
end 

#FVA: this double definition is so that MaxMinPlus can use it. 
"""
    ⊗̇(x::EntropyDoubleCSemifield{T,τ}, y::EntropyDoubleCSemifield{T,τ}) where {T,τ}

Upper multiplication in double complete semifields, including results for the top which is absorptive.
"""
⊗̇(x::EntropyDoubleCSemifield{T,τ}, y::EntropyDoubleCSemifield{T,τ}) where {T,τ} = oplusu(x, y)


""" 
    struct MaxMinPlusSemifield{T} <: DoubleCSemifield{T}

The double complete semifield that is aligned with max-plus.
"""
const MaxMinPlusSemifield{T} = EntropyDoubleCSemifield{T,Inf} where T

"""
    ⊕(x::S, y::S) where S<:MaxMinPlusSemifield = S(max(val(x), val(y)))
""" 
⊕̧(x::S, y::S) where S<:MaxMinPlusSemifield = S(max(val(x), val(y)))

"""
    ⊕̇(x::S, y::S) where S<:MaxMinPlusSemifield = S(min(val(x), val(y)))
"""
⊕̇(x::S, y::S) where S<:MaxMinPlusSemifield = S(min(val(x), val(y)))

#FVA: for some reason the inference over types is incomplete and 
#=
"""
    const TropicalDoubleCSemifield{T} = EntropyDoubleCSemifield{T,-Inf} where T

Completed tropical or min-plus semifield: ``(\\mathbb{R} \\cup \\{- \\infty \\}, min, +, \\infty, 0, -\\infty)``.

This is the completed entropy semifield with ``\\tau = -\\infty``.
"""
const TropicalDoubleCSemifield{T} = EntropyDoubleCSemifield{T,-Inf} where T
const MinPlusDoubleCSemifield{T} = TropicalDoubleCSemifield{T} where T

⊕(x::S, y::S) where S<:TropicalDoubleCSemifield = S(min(val(x), val(y)))
#⊗(x::S, y::S) where S<:TropicalDoubleCSemifield = S(val(x) + val(y))#Not changed from the original
Base.zero(S::Type{<:TropicalDoubleCSemifield{T}}) where T = S(Inf)#Shorter than the original
#Base.one(S::Type{<:TropicalDoubleCSemifield{T}}) where T = S(zero(T))#Not changed from the original
top(S::Type{<:TropicalDoubleCSemifield{T}}) where T = S(-Inf)#Shorter than the original
#Base.inv(x::S) where S<:TropicalDoubleCSemifield = S(-val(x))
#∂sum(z::S, x::S) where S<:TropicalDoubleCSemifield = valtype(S)(x == z)
#∂rmul(x::S, a::S) where S<:TropicalDoubleCSemifield = valtype(S)(1)#FVA: not sure of this
#∂lmul(a::S, x::S) where S<:TropicalDoubleCSemifield = valtype(S)(1)#FVA: not sure of this

"""
    const ArcticDoubleCSemifield{T} = EntropyDoubleCSemifield{T,Inf} where T

Complete arctic or max-plus semifield: ``\\langle \\mathbb{R} \\cup \\{\\infty \\}, max, +, -\\infty, 0, +\\infty \\rangle``.
This is the complete entropy semifield with ``\\tau = \\infty``.

"""
const ArcticDoubleCSemifield{T} = EntropyDoubleCSemifield{T,Inf} where T#Technically correct, but a mouthful
const MaxPlusDoubleCSemifield{T} = ArcticDoubleCSemifield{T} where T#Preferred named
⊕(x::S, y::S) where S<:ArcticDoubleCSemifield = S(max(val(x), val(y)))
#⊗(x::S, y::S) where S<:ArcticDoubleCSemifield = S(val(x) + val(y))#Not changed from the original
Base.zero(S::Type{<:ArcticDoubleCSemifield{T}}) where T = S(-Inf)#Shorter than the original
#Base.one(S::Type{<:ArcticDoubleCSemifield{T}}) where T = S(zero(T))#Not changed from the original
top(S::Type{<:ArcticDoubleCSemifield{T}}) where T = S(Inf)#Shorter than the original
#Base.inv(x::S) where S<:ArcticDoubleCSemifield = S(-val(x))#Not changed from the original
#∂sum(z::S, x::S) where S<:ArcticDoubleCSemifield = valtype(S)(x == z)#FVA:  WRONG after the original author of Semirings. 
#∂rmul(x::S, a::S) where S<:ArcticDoubleCSemifield = valtype(S)(1)#FVA: not sure of this
#∂lmul(a::S, x::S) where S<:ArcticDoubleCSemifield = valtype(S)(1)#FVA: not sure of this

=#
end#DoubleCSemifields
