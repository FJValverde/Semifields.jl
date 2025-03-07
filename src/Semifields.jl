# SPDX-License-Identifier: CECILL-2.1
"""
        Semifields

A module to capture the concept of a semifield and provide its most notorious examples.
"""
module Semifields

using ChainRulesCore
import LogExpFunctions: logaddexp

export
    # Main API.
    ⊕, #\oplus
    ⊗, #\otimes
    #inv, #conjugation
    top, 
    val,

    # Partial derivative functions.
    ∂sum,
    ∂rmul,
    ∂lmul,

    # Concrete types.
    Semifield,
    ArcticSemifield,
    BoolSemifield,
    LogSemifield,
    # ProbSemifield,
    TropicalSemifield


###############################################################################
# Abstract semifield type.
include("main_API.jl")


###############################################################################
# Concrete semifield types.
include("concrete_semifields.jl")

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

function ChainRulesCore.rrule(::typeof(val), a::Semifield)
    b = val(a)
    function val_pullback(b̄)
        ZeroTangent(), Tangent{typeof(a)}(; val=b̄)
    end
    b, val_pullback
end


function ChainRulesCore.rrule(S::Type{<:Semifield}, a)
    b = S(a)
    function semifield_pullback(b̄)
        ZeroTangent(), b̄.val
    end
    b, semifield_pullback
end

end#module Semifields
