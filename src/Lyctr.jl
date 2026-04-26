module Lyctr

using Oscar
import Oscar: evaluate, polynomial_ring
using FusionRings
using PicoSAT

include("general_functions.jl")
include("Gideon/Gideon.jl")
include("Harrow/Harrow.jl")

function __init__()
    # GLOBAL VARIABLES
    # Abreviations for commonly used fields
end

end
