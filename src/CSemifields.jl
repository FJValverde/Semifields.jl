using Reexport
@reexport module CompleteSemifields
"""
    CompleteSemifields -- Complete Semifields

A module to define complete semifields as pairs of two semifields
related by inversion.
"""

export 
    top, ⊤,
    ⊥, #an alternate name for the bottom element. 

    CompleteSemifield,#the basic abstract type.
    TernaryCSemifield#The basic semifield included in all the rest.
#=
    EntropyCSemifield,#The pair of semifields of logarithmic ranges.
    MaxMinPlus,#The pair of semifields of tropical ranges.
    #NaturalCSemifield,#The pair of semifields of natural ranges.
    MaxMinTimes
=#
# CAVEAT! The complete Semifields are never subtypes of the 
# concrete semifields

#using ChainRulesCore
#import LogExpFunctions: logaddexp
#using Reexport
#@reexport using ..Semifields#This should not be needed, since Semirings is re-exported by Semifields.

import ..Semirings: ⊕,⊗,val, valtype,∂sum,∂rmul,∂lmul, Semiring
import ..GenericSemifields: Semifield

"""
    CompleteSemifield

A type for abstract, generic complete semifields.
- Tp is the type of the semifield order-aligned with this semifield.
- Tn is the type of the dually order-aligned semifield. 
"""
abstract type CompleteSemifield{T} <: Semifield{T}
end

"""
    TernarySemifield <: CompleteSemifield{Some(Bool)}
"""
struct TernarySemifield <: CompleteSemifield{Some{Bool}}
    val::Some{Bool}
end
#Base.zero(::Some{Bool}) = TernarySemifield(false)
Base.zero(::Type{<:TernarySemifield}) = TernarySemifield(Some(false))
top(::Type{<:TernarySemifield}) = TernarySemifield(Some(true))
Base.one(::Type{<:TernarySemifield}) = TernarySemifield(nothing)
Base.inv(x::TernarySemifield) = isnothing(x.val) ? x : TernarySemifield(!x.val)


end #module CompleteSemifields