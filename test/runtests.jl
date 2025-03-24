
using ChainRulesCore
using LogExpFunctions
using Semifields
using Test

begin
include("semiringtests.jl")
end
begin
include("semifieldtests.jl")
end
begin
include("csemifieldtests.jl")
end
begin
include("doublecsemifieldtests.jl")
end     

#include("matrixtests.jl")

