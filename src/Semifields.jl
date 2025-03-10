# SPDX-License-Identifier: CECILL-2.1
"""
        Semifields

A module to capture the concept of a semifield and provide its most notorious examples.
"""
module Semifields

#FVA all main imports to locate them. 
using Reexport
using ChainRulesCore
import LogExpFunctions: logaddexp

export
    Semirings,#The base module on semirings.
    GenericSemifields,#The module on generic semifields.
    CompleteSemifields#The module on complete semifields.

#FVA the following module is completely re-exported. 
include("semirings.jl")#An in-house re-write of L.Ondel Yang's code. 

#using Reexport
#@reexport 
module GenericSemifields#FVA: cannot be called Semimodules due to parent
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
    # ProbSemifield,#Aka Energy semiring, natural semiring, etc.
    TropicalSemifield

###############################################################################
# Abstract semifield type.
include("main_API.jl")
###############################################################################
# Concrete semifield types.
include("concrete_semifields.jl")
end #module GenericSemifields
#
#FVA: the following shamelessly copies the matrixsemiring.jl file by L. Ondel.
#module MatrixSemirings

#export MatrixSemiring, +, *
#include("matrixsemiring.jl")#may be not needed!

#end
#using Reexport
#@reexport using .Semirings#defined in semirings.jl

#Reexport.@reexport MatrixSemirings
#export MatrixSemifield, +, *, inv, top
#include("matrixsemifield.jl")

include("CSemifields.jl")

end#module Semifields
