module Lyctr

using Oscar
using FusionRings
using PicoSAT

include("general_functions.jl")
include("Gideon/Gideon.jl")
include("Harrow/Harrow.jl")

function __init__()
    # GLOBAL VARIABLES
    # Abreviations for commonly used fields
    ℚ = QQ;
    ℚb       = algebraic_closure(ℚ)
    ℚab, ζ   = abelian_closure(ℚ; sparse = false)
    sℚab, ζs = abelian_closure(ℚ; sparse = true)
end

end
