module Lyctr

using Oscar
import Oscar: evaluate, polynomial_ring, sparse_matrix, snf
using FusionRings
import FusionRings: fusion_ring, multiplication_table, frl, fawc
using PicoSAT
using TensorCategories
import TensorCategories: F_symbols, R_symbols, P_symbols, anyonwiki, multiplication_table

include("general_functions.jl")
include("Gideon/Gideon.jl")
include("Harrow/Harrow.jl")

function __init__()
    # GLOBAL VARIABLES
    # Abreviations for commonly used fields
end

end
