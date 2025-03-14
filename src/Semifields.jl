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

#FVA the following module is completely re-exported. 
include("generic_semifields.jl")#A module to define generic semifields.

# The following module is completely re-exported.
include("CSemifields.jl")#A module to define complete semifields.

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


end#module Semifields
