"""
    CompleteSemifield

A type for abstract, generic complete semifields.
- Tp is the type of the semifield order-aligned with this semifield.
- Tn is the type of the dually order-aligned semifield. 
"""
abstract type CompleteSemifield{Tp,Tn} <: Semifield{Tp}
end

"""
    TernarySemifield <: CompleteSemifield{Some(Bool)}
"""
struct TernarySemifield <: CompleteSemifield{Some(Bool)}
    val::Some(Bool)
end
Base.zero(Some(Bool)) = TernarySemifield(false)
Base.zero(::Type{<:TernarySemifield}) = TernarySemifield(Some(false))
top(::Type{<:TernarySemifield}) = TernarySemifield(Some(true))
Base.one(::Type{<:TernarySemifield}) = TernarySemifield(nothing)
"""
    TernaryCSemifield{T} = CompleteSemifield{T,T} where T
"""
type TernaryCSemifield{T} = CompleteSemifield{BooleanSemifield,BooleanSemifield}
end

"""
    EntropyCSemifield{T,R} = CompleteSemifield{LogSemifield{T,R}, LogSemifield{T,-R}} where {T,R}
"""
const EntropyCSemifield{T,R} = 
    CompleteSemifield{LogSemifield{T,R}, LogSemifield{T,-R}} where {T,R}

const MaxMinPlus{T} = 
    CompleteSemifield{TropicalSemifield{T}, ArcticSemifield{T}} where {T,R}

const MaxMinTimes{T} = 
    CompleteSemifield{MaxTimesSemifield{T}, MinTimesSemifield{T}} where {T,R}