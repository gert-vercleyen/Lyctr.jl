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

export QQb, QQab, ζ, QQabs, QQemb 
function __init__()
  global QQb      = algebraic_closure(QQ)
  global QQab, ζ  = abelian_closure(QQ)
  # Following is necessary to work with trivial embeddings from ℚ → ℂ
  global QQabs, _ = rationals_as_number_field()
  global QQemb    = hom( QQabs, QQb, QQb(1) )

  # GLOBAL VARIABLES
  # Abreviations for commonly used fields

  # importing data 
  global datadir = joinpath(@__DIR__,"data")
  fns = readdir(datadir)

  if "split_number_data" ∉ fns
    println("Dataset of numbers not yet optimized. Optimizing for future use.")

    # Create directory for numbers
    splitdatapath = joinpath(datadir, "split_number_data/")
    qqbpath       = joinpath(splitdatapath,"algebraic_numbers")
    asnpath       = joinpath(splitdatapath,"abstract_numbers")

    mkdir(splitdatapath)
    mkdir(qqbpath)
    mkdir(asnpath)

    println("Importing numbers.")
    # Import qqb numbers from big files
    qqb_nums = begin
      ids  = Oscar.load(joinpath(datadir, "qqb_ids.mrdi"))
      nums = Oscar.load(joinpath(datadir, "qqb_vals.mrdi"))

      [ids[i], nums[i] for i in 1:length(ids)]
    end
    #=
    absnumfieldelems = begin
      nothing
    end 

    embeddings = begin
      nothing
    end
    =#
    println("Exporting numbers separately.")
    # Export qqb numbers
    function exportnum(tuple)
      dir = joinpath(datadir, "split_number_data/", tuple[1]*".mrdi")
      return Oscar.save(path, tuple[2])
    end

    exportnum.(qqb_nums)
    println("Dataset is optimized.")
  end

  # start an empty qqb_dict which grows on demand 
  global qqb_dict = Dict{String, QQBarFieldElem}()
end

end
