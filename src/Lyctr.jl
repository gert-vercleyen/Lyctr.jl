module Lyctr

using  Oscar
import Oscar: evaluate, polynomial_ring, sparse_matrix, snf, ring, embedding, base_field

using  FusionRings
import FusionRings: fusion_ring, multiplication_table, frl, fawc

using  PicoSAT

using  TensorCategories
import TensorCategories:
  F_symbols, R_symbols, P_symbols, anyonwiki, multiplication_table

import Base.values, Base.merge, Base.length

using ProgressMeter

using JSON

include("general_functions.jl")
include("Gideon/Gideon.jl")
include("Harrow/Harrow.jl")

function __init__()
  global QQb     = algebraic_closure(QQ)
  global QQab, ζ = abelian_closure(QQ)
    # GLOBAL VARIABLES
    # Abreviations for commonly used fields
end

end
