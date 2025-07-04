# SPDX-License-Identifier: CECILL-2.1
"""
        Semifields

A module to capture the concept of a semifield and provide its most notorious example.

Semifields are semirings whose multiplication is invertible, so we also provide semirings
as a submodule.

Complete semifields have a top element, and are defined as pairs of two semifields related 
by inversion.
"""
module Semifields

#FVA all main imports to locate them. 
using Reexport
using ChainRulesCore
import LogExpFunctions: logaddexp

export
    Semirings,#The base module on semirings. Re-exported.
    MatrixSemirings,#The module on matrix semirings. Re-exported.
    #GenericSemifields,#The module on generic semifields.
    CSemifields,#The module on complete semifields. Not re-exported.
    DoubleCSemifields#The module on double complete semifields. Re-exported.

#FVA the following module is completely re-exported. 
include("semirings.jl")#An in-house re-write of L.Ondel Yang's code. 

#FVA the following module is completely re-exported. 
include("generic_semifields.jl")#A module to define generic semifields.

# The following module is completely re-exported.
#using Reexport
#@reexport using CSemifields
include("CSemifields.jl")#A module to define complete semifields.

#FVA: use the CSemifields module to define the double complete semifields.
include("DoubleCSemifields.jl")#A module to define complete semifields.

#FVA: the following shamelessly copies the matrixsemiring.jl file by L. Ondel.
#module MatrixSemirings

#export MatrixSemiring, +, *
#FVA: the following module is completely re-exported.
# It is a test on Julia's capabilities to define a semiring on matrices.
# Note that this is not a complete semiring, but a semiring on matrices.
# Furthermore, this models left of right actions on matrix semimodules, TBD
include("matrixsemiring.jl")

end#module Semifields
