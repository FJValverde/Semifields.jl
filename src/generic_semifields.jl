using Reexport
@reexport module GenericSemifields#FVA: cannot be called Semimodules due to parent
"""
    GenericSemifields

A module to define generic and concrete semifields. 
"""

export
    #⊕, #\oplus#FVA: already defined on semirings.jl
    #⊗, #\otimes#FVA: already defined on semirings.jl
    inv, #conjugation
    top, #this belongs to complete semifields. 
    #val,#FVA: already defined on semirings.jl
    #valtype,#FVA: already defined on semirings.jl

    # Partial derivative functions.
    #∂sum,#FVA: already defined on semirings.jl
    #∂rmul,#FVA: already defined on semirings.jl
    #∂lmul,#FVA: already defined on semirings.jl

    # Concrete types.
    Semifield,
    BoolSemifield,
    EntropySemifield,
    ArcticSemifield,
    TropicalSemifield
    # ProbSemifield,#Aka Energy semiring, natural semiring, etc.

###############################################################################
# Abstract semifield type.
include("main_API.jl")
###############################################################################
# Concrete semifield types.
include("concrete_semifields.jl")
end #module GenericSemifields
