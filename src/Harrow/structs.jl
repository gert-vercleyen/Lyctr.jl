
#┏━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━┓
#┃                                      FL                                         ┃
#┗━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━┛

# Shorthand for F-label
FL = NTuple{10, Int64}

#┏━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━┓
#┃                                      RL                                         ┃
#┗━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━┛

# Shorthand for R-label
RL = NTuple{5, Int64}

#┏━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━┓
#┃                                      PL                                         ┃
#┗━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━┛

# Shorthand for P-label
PL = NTuple{1, Int64}

#┏━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━┓
#┃                                      GL                                         ┃
#┗━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━┛

# Shorthand for G-label
GL = NTuple{5, Int64}

#┏━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━┓
#┃                                 SkelFusCatSymb                                  ┃
#┗━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━┛

abstract type AbstractSkelFusCatSymb end

ASFCS = AbstractSkelFusCatSymb

struct SkelFusCatSymb <: AbstractSkelFusCatSymb
  # dict mapping indices to values
  dict

  # polynomial ring in symbols. We store this since, once all values of th
  # symbols are known, there is no way to retrieve this ring via the
  # parent function
  polyring::Ring

  # dictionary mapping generators of polyring to labels
  # This keeps track of the relation between labels and generators
  # so that fixing values of dict does not destroy the relationship
    
  idict

  # the base_field is available via the polyring but if we want the field
  # to be embedded we need to store this separately.
  # Note that, when working over QQb or QQab, no embedding is
  # needed as these fields are embedded by default
  embedding::CEmbOrNothing

  # Only store multiplication table of fusion ring since
  # storing the whole ring takes space.
  mt::Array{Int64, 3}

  #inminfield is true if symbols are expressed over a minimal field
  inminfield::Union{Bool, Missing}

  #inunitarygauge is true if the symbols are expressed in a unitary gauge
  inunitarygauge::Union{Bool, Missing}

end

SFCS = SkelFusCatSymb

function skel_fus_cat_symb(
  d::Dict,
  idict,
  polring::Ring,
  embedding::CEmbOrNothing,
  mt;
  inminfield = missing,
  inunitarygauge = missing,
)
  return SkelFusCatSymb(d, polring, idict, embedding, mt, inminfield, inunitarygauge)
end

function skel_fus_cat_symb(
  d::Dict,
  polring::Ring,
  idict,
  mt;
  inminfield = missing,
  inunitarygauge = missing
  )

  F = base_ring(polring)
  emb = F == QQ ? QQemb : nothing
  
  return SkelFusCatSymb(d, polring, idict, emb, mt, inminfield, inunitarygauge)
  
end

    

export dict

"""dict(fs::SkelFusCatSymb) returns a dictionary mapping the labels of the symbols to their values
"""
function dict(fs::ASFCS)
  return fs.dict
end

export labels

"""labels(fs::SkelFusCatSymb) returns the labels of the symbols in lexicographic order
"""
function labels(fs::ASFCS)
  return (sort∘collect∘keys∘dict)(fs)
end

export values

"""values(fs::SkelFusCatSymb) returns the values of the symbols, in lexicographic order on the labels
"""
function values(fs::ASFCS)
  d = dict(fs)
  labls = labels(fs)
  return [ d[l] for l in labls ] 
end

export invdict

function invdict(fs::ASFCS) 
  return fs.idict
end 

export value

function value(fs::ASFCS)::Function
  return x -> dict(fs)[x]
end

function value(fs::ASFCS, x)::RingElem
  return dict(fs)[x]
end

function value(fs::ASFCS, v::AbstractVector)::Vector{RingElem}
  return [ dict(fs)[x] for x in v ]
end

function value(fs::ASFCS, m::AbstractMatrix)::Matrix{RingElem}
  return [ dict(fs)[m[i,j]] for i in 1:nrows(m), j in 1:ncols(m) ]
end

export polynomial_ring

function polynomial_ring(fs::ASFCS)
  return fs.polyring
end

export base_field

function base_field(fs::ASFCS)
  return base_ring(polynomial_ring(fs))
end

export embedding

function embedding(fs::ASFCS)
  return fs.embedding
end
  
export multiplication_table

function multiplication_table(fs::ASFCS)
  return fs.mt
end

export fusion_ring

"""fusion_ring(s::SFCS) returns the fusion ring belonging to
the fusion category with symbols s
"""
function fusion_ring(fs::ASFCS)
  return fusion_ring(multiplication_table(fs))
end

export inminfield

function inminfield(fs::ASFCS)
  return fs.inminfield
end

export inunitarygauge

function inunitarygauge(fs::ASFCS)
  return fs.inunitarygauge
end

function length(fs::ASFCS)
  return length(dict(fs))
end


# allvaluesknown is true if all_values of the dict of symbols belong to
# base_ring(polyring)
export allvaluesknown

function allvaluesknown(fs::ASFCS)
  return all( is_constant, values(fs) )
end

function evaluate( p::RingElem, fs::ASFCS )
  is_constant(p) && return p

  varind = var_indices(p)
  R      = parent(p)
  vls    = [ dict(fs)[invdict(fs)[R[i]]] for i in varind ]

  
  #@info "polynomial: " p
  #@info "parent of polynomial: " parent(p)
  #@info "base field polynomial: " base_ring(p)
  #@info "indices of pol :" varind 
  #@info "values for evaluation: " vls
  #@info "parents of values: " parent.(vls)
  #@info "base fields of values: " base_ring.(vls)
  
  
  evaluate(p,varind,vls) 
end
  
function evaluate(vp::Vector{RingElem},fs::ASFCS)
  [ evaluate( p, fs ) for p in vp ]
end

export merge 

"""merge(symb1::ASFCS,symb2::ASFCS)::ASFCS

returns a ASFCS struct whose polynomial ring contains elements of
 of symb1 and symb2. Only works atm if base_field(symb1) == base_field(symb2)
"""

function merge(symb1::ASFCS,symb2::ASFCS)
  base_field(symb1) != base_field(symb2) && error("Number fields of symbols must be the same")
  embedding(symb1) != embedding(symb2) && error("Embeddings of symbols must be the same")
  mtab = multiplication_table
  mtab(symb1) != mtab(symb2) && error("Multiplication tabels of fusion rings of symbols must be the same")

  R1 = polynomial_ring(symb1)
  R2 = polynomial_ring(symb2)
  R, vrs, i1, i2 = merge( R1, R2 )

  # embed dictionaries
  ed1 = Dict( (k, i1(v) ) for (k,v) ∈ dict(symb1) )
  ed2 = Dict( (k, i2(v) ) for (k,v) ∈ dict(symb2) )

  newdict = Base.merge(ed1,ed2)

  idict = merge(invdict(symb1),invdict(symb2))
  
  emb   = embedding(symb1)

  mt    = mtab(symb1)

  imf1 = inminfield(symb1)
  imf2 = inminfield(symb2)
  imf  = ismissing(imf1) || ismissing(imf2) ? missing : imf1 && imf2

  iug1 = inunitarygauge(symb1)
  iug2 = inunitarygauge(symb2)
  iug  = ismissing(iug1) || ismissing(iug2) ? missing : iug1 && iug2
  
  return skel_fus_cat_symb(newdict,idict,R,emb,mt;inminfield=imf,inunitarygauge=iug)
end

function Base.show(io::IO, s::ASFCS )
  d = dict(s)
  l = length(d)
  sortedsymb = sort(collect(d), by = (x->x[1]))

  headerstring = 
    "Collection of " * 
    string(l) * 
    " labeled skeletal fusion category symbols over " *
    string(base_field(s))
  
  embstring =
    if embedding(s) === nothing
      "\n"
    else
      ", embedded via\n" * string(embedding(s)) * "\n"
    end
  symbstring =
    if l <= 13
      "Symbols:\n" * join( string.(sortedsymb) , "\n" )
    else
"\nSymbols:\n" *
      join( string.(sortedsymb)[1:6], "\n" )*
        "\n⋮\n" *
        join( string.(sortedsymb)[end-5:end], "\n" )
    end
  print(io, headerstring * embstring * symbstring)
      
end


#┏━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━┓
#┃                                    FSymbols                                     ┃
#┗━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━┛

abstract type AbstractFSymbols <: AbstractSkelFusCatSymb end

struct FSymbols <: AbstractFSymbols
  # dict mapping indices to values
  dict

  # polynomial ring in symbols
  polyring

  idict

  embedding::CEmbOrNothing

  # multiplication table of fusion ring. Storing the whole ring
  # takes space
  mt::Array{Int64, 3}

  #inminfield is true if symbols are expressed over a minimal field
  inminfield::Union{Bool, Missing}

  #inunitarygauge is true if the symbols are expressed in a unitary gauge
  inunitarygauge::Union{Bool, Missing}
end

function F_symbols(
  d::Dict{FL, T},
  polring::Ring,
  idict,
  embedding,
  mt;
  inminfield = missing,
  inunitarygauge = missing,
) where {T <: Union{ RingElem, FieldElem} }
  return FSymbols(d, polring, idict, embedding, mt, inminfield, inunitarygauge)
end

function F_symbols(
  d::Dict{FL, T},
  polring::Ring,
  idict,
  mt;
  inminfield = missing,
  inunitarygauge = missing
) where {T <: Union{ RingElem, FieldElem} }
  
  F = base_ring(polring)
  emb = F == QQ ? QQemb : nothing
  
  return F_symbols(d, polring, idict, emb, mt, inminfield, inunitarygauge )
end

#┏━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━┓
#┃                                 InverseFSymbols                                 ┃
#┗━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━┛

# InverseFSymbols has exactly the same property as FSymbols but a different
# name and printing. It only exists to avoid confusion.
abstract type AbstractInverseFSymbols <: AbstractSkelFusCatSymb end

struct InverseFSymbols <: AbstractInverseFSymbols
  # dict mapping indices to values
  dict

  # polynomial ring in the variables F[1], ..., F[n]
  polyring

  idict

  embedding::CEmbOrNothing

  # multiplication table of fusion ring. Storing the whole ring
  # takes space
  mt::Array{Int64, 3}

  #inminfield is true if F-symbols are expressed over a minimal field
  inminfield::Union{Bool, Missing}

  #inunitarygauge is true if the F-symbols are expressed in a unitary gauge
  inunitarygauge::Union{Bool, Missing}
end

function inverse_F_symbols(
  d::Dict{FL, T},
  R::Ring,
  idict,
  embedding::CEmbOrNothing, 
  mt;
  inminfield = missing,
  inunitarygauge = missing,
) where {T <: Union{RingElem,FieldElem} }
  return InverseFSymbols(d, R, idict, embedding, mt, inminfield, inunitarygauge )
end

function inverse_F_symbols(
  d::Dict{FL, T},
  R::Ring,
  idict,
  mt;
  inminfield = missing,
  inunitarygauge = missing,
) where {T <: Union{ RingElem, FieldElem} }
  
  F = base_ring(R)
  emb = F == QQ ? QQemb : nothing
  
  return inverse_F_symbols(d, R, idict, emb, mt, inminfield, inunitarygauge )
end

#┏━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━┓
#┃                                    RSymbols                                     ┃
#┗━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━┛

abstract type AbstractRSymbols <: AbstractSkelFusCatSymb end

struct RSymbols <: AbstractRSymbols
  # dict mapping indices to values
  dict

  # polynomial ring in symbols
  polyring::Ring

  idict

  embedding::CEmbOrNothing

  # multiplication table of fusion ring. Storing the whole ring
  # takes space
  mt::Array{Int64, 3}

  #inminfield is true if symbols are expressed over a minimal field
  inminfield::Union{Bool, Missing}

  #inunitarygauge is true if the symbols are expressed in a unitary gauge
  inunitarygauge::Union{Bool, Missing}
end

function R_symbols(
  d::Dict{RL, T},
  R::Ring,
  idict,
  embedding,
  mt;
  inminfield = missing,
  inunitarygauge = missing
) where {T <: RingElem}
  return RSymbols(d, R, idict, embedding, mt, inminfield, inunitarygauge)
end

function R_symbols(
  d::Dict{RL, T},
  R::Ring,
  idict,
  mt;
  inminfield = missing,
  inunitarygauge = missing,
) where {T <: Union{ RingElem, FieldElem} }
  
  F = base_ring(R)
  emb = F == QQ ? QQEmb : nothing
  
  return R_symbols(d, R, idict, emb, mt, inminfield, inunitarygauge )
end

#┏━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━┓
#┃                                    InverseRSymbols                                     ┃
#┗━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━┛

abstract type AbstractInverseRSymbols <: AbstractSkelFusCatSymb end

struct InverseRSymbols <: AbstractInverseRSymbols
  # dict mapping indices to values
  dict

  # polynomial ring in symbols
  polyring::Ring

  idict

  embedding::CEmbOrNothing
  # multiplication table of fusion ring. Storing the whole ring
  # takes space
  mt::Array{Int64, 3}

  #inminfield is true if symbols are expressed over a minimal field
  inminfield::Union{Bool, Missing}

  #inunitarygauge is true if the symbols are expressed in a unitary gauge
  inunitarygauge::Union{Bool, Missing}
end

function inverse_R_symbols(
  d::Dict{RL, T},
  R::Ring,
  idict,
  embedding,
  mt;
  inminfield = missing,
  inunitarygauge = missing,
  allvaluesknown = missing,
) where {T <: RingElem}
  return InverseRSymbols(d, R, idict, embedding, mt, inminfield, inunitarygauge)
end

function inverse_R_symbols(
  d::Dict{FL, T},
  R::Ring,
  idict,
  mt;
  inminfield = missing,
  inunitarygauge = missing,
) where {T <: Union{ RingElem, FieldElem} }
  
  F = base_ring(R)
  emb = F == QQ ? QQEmb : nothing
  
  return inverse_R_symbols(d, R, idict, emb, mt, inminfield, inunitarygauge )
end
    
#┏━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━┓
#┃                                    PSymbols                                     ┃
#┗━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━┛

abstract type AbstractPSymbols <: AbstractSkelFusCatSymb end

struct PSymbols <: AbstractPSymbols
  # dict mapping indices to values
  dict

  # polynomial ring in symbols
  polyring

  idict

  embedding::CEmbOrNothing

  # multiplication table of fusion ring. Storing the whole ring
  # takes space
  mt::Array{Int64, 3}

  #inminfield is true if symbols are expressed over a minimal field
  inminfield::Union{Bool, Missing}

  #inunitarygauge is true if the symbols are expressed in a unitary gauge
  inunitarygauge::Union{Bool, Missing}
end

function P_symbols(
  d::Dict{PL, T},
  R::Ring,
  idict,
  embedding::CEmbOrNothing,
  mt;
  inminfield = missing,
  inunitarygauge = missing
) where {T <: RingElem}
  return PSymbols(d, R, idict, embedding, mt, inminfield, inunitarygauge)
end

function P_symbols(
  d::Dict{PL, T},
  R::Ring,
  idict,
  mt;
  inminfield = missing,
  inunitarygauge = missing,
) where {T <: Union{ RingElem, FieldElem} }
  
  F = base_ring(R)
  emb = F == QQ ? QQEmb : nothing
  
  return P_symbols(d, R, idict, emb, mt, inminfield, inunitarygauge )
end

#┏━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━┓
#┃                                    GSymbols                                     ┃
#┗━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━┛

abstract type AbstractGSymbols <: AbstractSkelFusCatSymb end

struct GSymbols <: AbstractGSymbols
  # dict mapping indices to values
  dict

  # polynomial ring in symbols
  polyring

  idict

  embedding::CEmbOrNothing

  # multiplication table of fusion ring. Storing the whole ring
  # takes space
  mt::Array{Int64, 3}
  
  #inminfield is true if symbols are expressed over a minimal field
  inminfield::Union{Bool, Missing}

  #inunitarygauge is true if the symbols are expressed in a unitary gauge
  inunitarygauge::Union{Bool, Missing}
end

function G_symbols(
  d::Dict{RL, T},
  R::Ring,
  idict,
  embedding::CEmbOrNothing,
  mt;
  inminfield = missing,
  inunitarygauge = missing
) where {T <: RingElem}
  return GSymbols(d, R, idict, embedding, mt,  inminfield, inunitarygauge)
end


function G_symbols(
  d::Dict{FL, T},
  R::Ring,
  idict,
  mt;
  inminfield = missing,
  inunitarygauge = missing,
) where {T <: Union{ RingElem, FieldElem} }
  
  F = base_ring(R)
  emb = F == QQ ? QQEmb : nothing
  
  return G_symbols(d, R, idict, emb, mt, inminfield, inunitarygauge )
end

#┏━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━┓
#┃                                   Symmetries                                    ┃
#┗━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━┛


abstract type AbstractFusCatSymbSymmetries end

# TODO: finish implementation
struct FusCatSymbSymmetries <: AbstractFusCatSymbSymmetries
  # symbols on which symmetries work
  fuscatsymb::SkelFusCatSymb

  # fraction field of gauge-symbols and symbols on which gauge-transforms work
  fracfield 

  # dictionary that maps symbols to their transformed form
  action::Dict{FieldElem,FieldElem}
end

function fus_cat_symb_symmetries( symb::SkelFusCatSymb, ff, dict )
  FusCatSymbSymmetries( symb, ff, dict )
end


#┏━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━┓
#┃                                SkeletalFusionCat                                ┃
#┗━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━┛

abstract type AbstractSkeletalFusionCat end

export SkeletalFusionCat

struct SkeletalFusionCat <: AbstractSkeletalFusionCat
  # Fusion ring identifier
  fr::FusionRing

  # ID
  id::String

  #uuid

  # F-symbols, R-symbols, P-symbols (a.k.a. pivotal coeffs)
  fsymb::FSymbols
  rsymb::Union{RSymbols, Missing, Nothing}
  psymb::Union{PSymbols, Missing, Nothing}
  # G-symbols: symbols for gauge transforms
  gsymb::GSymbols

  # field over which F-, R-, and P-symbols are defined
  basefield::Field

  # image of embedding of generator of abstract number field in field
  # of algebraic numbers: QQBar.
  # nothing if field is QQb, QQab, finite, p-adic, or abstract
  embedding::Union{QQBarFieldElem,Nothing}

  # inminfield is true if ALL symbols are expressed over a minimal field
  inminfield::Union{Bool, Missing}

  # minimalfields is a dict that maps the strings
  # "F", "FR", "FP", and "FPR" to the mimimal fields for the F-symbols
  # F- and R-symbols, F- and P-symbols and F-, P-, and R-symbols
  minimalfields::Union{Dict{String,Field},Missing}

  # Inunitarygauge is
  # * true if ALL symbols are expressed in a unitary gauge
  # * missing if info is not available
  # * nothing if unitarity doesn't make sense due to field of symbols
  inunitarygauge::Union{Bool, Missing}

  # names of the fusion category (regular names with an a- in front or so)
  nms::Union{Missing, Dict{String,Union{Vector{String},Missing}}}

  #tex names of the cat
  texnms::Union{Missing, Dict{String,Union{Vector{String},Missing}}}

  # non-trivial sub fusion cats
  # sub fusion cats are given in a dict that maps identifiers
  # to embeddings
  ntsfc::Union{Missing,Dict{String,Vector{Int}}}

  # inverse gauge-split transform: sparse matrix that transforms
  # gauge invariants back to F, R and P symbols
  igst::Union{Missing,SMat}

  # gauge-split basis: tuple (I,D) of lists of words in formal
  # F, R, and P symbols where all el of I are gauge-invariant,
  # all el of D are gauge-dependent and any of the F, R, P symbols
  # can be uniquele written as a product of the words
  gsb::Union{Missing,Tuple{Vector{Any},Vector{Any}}}

  #realizations via some known constructions
  rlztns::Union{Missing,Dict{String,Any}}

  #is pivotal
  isp::Union{Missing,Bool}
  #is unitary
  isu::Union{Missing,Bool}
  #is braided
  isb::Union{Missing,Bool}
  #is modular
  ism::Union{Missing,Bool}
end

export skeletal_fusion_cat

"""skeletal_fusion_cat( fr, id, fs, abstractfield ) returns a fusion category over an abstract field
"""
function skeletal_fusion_cat(
  fusionring,
  id,
  fsymbols,
  basefield,
  embedding;
  rsymbols                   = missing,
  psymbols                   = missing,
  inminimalfield             = missing,
  minimalfields              = missing,
  inunitarygauge             = missing,
  names                      = missing,
  texnames                   = missing,
  nontrivialsubcategories    = missing,
  inversegaugesplittransform = missing,
  gaugesplitbasis            = missing,
  realizations               = missing,
  ispivotal                  = missing,
  isunitary                  = missing,
  isbraided                  = missing,
  ismodular                  = missing
  )

  nms = if ismissing(names)
    Dict(
      "quantum_group_like" => missing,
      "group_like" => missing,
      "physics" => missing,
      "miscelaneous" => missing
    )
  else
    names
  end
  
  texnms = if ismissing(texnames)
    Dict(
      "quantum_group_like" => missing,
      "group_like" => missing,
      "physics" => missing,
      "miscelaneous" => missing
    )
  else
    texnames
  end
  
  return SkeletalFusionCat(
    fusionring,
    id,
    fsymbols,
    rsymbols,
    psymbols,
    G_symbols(fusionring,field=basefield),
    basefield,
    embedding,
    inminimalfield,
    minimalfields,
    inunitarygauge,
    nms,
    texnms,
    nontrivialsubcategories,
    inversegaugesplittransform,
    gaugesplitbasis,
    realizations,
    ispivotal,
    isunitary,
    isbraided,
    ismodular
  )
end
