#┏━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━┓
#┃                          Exporting and importing numbers                        ┃
#┗━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━┛

# Exact numbers objects are quite heavy to load so we will load them all in
# a dictionary and provide functions to convert numbers to keys and vice versa

# path where data is stored


# In the latter, we often need to construct a string from a polynomial
function polstring(x::RingElem)
  return join( string.( collect( coefficients(x) ) ), "_" )
end

#┏━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━┓
#┃                     Exporting and importing QQBarFieldElems                     ┃
#┗━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━┛

# Generate unique ID for a QQBarFieldElem
export qqb_id

function qqb_id(x::QQBarFieldElem)
  mp = minimal_polynomial(x)

  plstring  = polstring(mp)
  numstring = string(rootnum(x))

  return plstring * "__" * numstring
end

function qqb_id(arr::Array{QQBarFieldElem})
  return qqb_id.(arr)
end

function rootnum(x::QQBarFieldElem)
  p   = minimal_polynomial(x)
  rts = roots(QQBar, p)
  sr  = sort(rts; by = root_sort_crit)
  return findfirst(y -> y == x, sr)
end

# This is the sort criterion for roots used by mathematica 
# and by the anyonwiki on 28/12/2025
function root_sort_crit(x)
  return (- Int(is_real(x)), real(x), imag(x))
end

function save_qqb_num(x::QQBarFieldElem)
  path = joinpath(@__DIR__, "data", "split_number_data", qqb_id(x)*".mrdi")
  return Oscar.save(path, x)
end

# Loading from library of qqb elements
# If number already loaded, use that one, otherwise load and
# add to qqb_dict

export from_qqb_id

function from_qqb_id(s::String)
  if haskey(qqb_dict, s)
    return qqb_dict[s]
  else
    fn = joinpath(datadir, "split_number_data", s*".mrdi")

    if isfile(fn)
      val = Oscar.load(fn)
      qqb_dict[s] = val
      return val
    else
      spl = split(s, "_")
      l   = length(spl)
      n   = parse(Int64, last(spl))
      cfs = ZZ.(parse.(Int64, spl[1:(end - 2)]))

      qqbval = sort(roots(QQBar, polynomial(ZZ, cfs)); by = root_sort_crit)[n]

      save_qqb_num(qqbval)

      return qqbval
    end
  end
end

function from_qqb_id(a::Vector)
  return from_qqb_id.(a)
end

function from_qqb_id(a::Array)
  return from_qqb_id.(a)
end

function from_qqb_id(a::Matrix)
  return from_qqb_id.(a)
end
export from_qqb_id


#┏━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━┓
#┃                     Exporting and importing complex floats                      ┃
#┗━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━┛


function reim(x::ComplexF64)::Vector{Float64}
  return [real(x), imag(x)]
end


#┏━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━┓
#┃                    Exporting and importing fusion categories                    ┃
#┗━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━┛


fusion_cat_info_text = """
    The following conventions are used in explaining the values of the dictionary:
    * While the data describes a (pivotal) (braided) fusion category, we will always refer to this category as a fusion category to save space.
    * A qqb_id is a string that uniquely describes an algebraic number. It is formatted as a list of n integers a0, ..., an, separated by underscores, followed by a double underscore, followed by an integer i: "a0_a1_..._an__i". Here a0 to an are the coefficients of a polynomial a0 + a1*x + ... + an*x^n and i denotes the i'th root of that polynomial. The indexing of roots of the polynomial takes the real roots first, in increasing order. Then come complex conjugate pairs of roots, sorted first by increasing real part and second by increasing complex part.
    * A af_id is a string that describes an element of an abstract field. It is represented as an array of n integers [a0, ..., an]. If x is the generator of the abstract field (stored under the key base_field as part of the data of a fusion category), then the vector represents the element a0 + a1*x + ... + an*x^n where n is the order of the field minus 1.
    * A ctf, or complex tuple of floats, is a representation of a complex floating point number a + i b by vector [ a, b ] where a and b are real floating point numbers.
    * If A is a JSON array then A[i] is its i'th element.
    * The mother category is the current category being represented by the JSON dictionary. When talking about subcategories, the mother category represents the category of which we're interested in its subcategories.
    * Every fusion category is uniquely identified by a uuid. By cat[cat_uuid] we mean the fusion category with uuid equal to cat_uuid
    * We work skeletally and therefore identify equivalence classes of objects with objects themselves. 
    * All stored fusion categories have a fixed order of their simple objects which equals the order of the elements of the fusion ring that is referenced. Therefore every simple object of a fusion category can be represented by a positive integer from 1 to r = rank(cat) we will often use elements and indices representing those elements interchangeably. By the i'th element of the cat[uuid] we mean the i'th element for the stored fusion category unless stated otherwise.
    * A value of null for any field means the data is missing unless states otherwise. In particular, it does not imply that the data doesn't exist.
    * A sof, or string of labels, is a string of integers seperated by underscores, i.e. of the form "a1_a2_..._an" where each ai is an integer. Its order is by definition the number of integers in the string.
    * By a P-symbol we mean a pivotal coefficient. (see "G. Vercleyen, On Low-Rank Multiplicity-Free Fusion Categories. National University of Ireland, Maynooth (Ireland), 2024." for the definition) 
    * The labels of the F-symbols, R-symbols, and P-symbols follow the conventions of "G. Vercleyen, On Low-Rank Multiplicity-Free Fusion Categories. National University of Ireland, Maynooth (Ireland), 2024." Where the labels of the symbol [F^{a,b,c}_d]^{(e,α,β)}_{(f,γ,δ)} are encoded as the string "a_b_c_d_e_α_β_f_γ_δ", the labels of [R^{a,b}_c]^{α}_{β} are encoded as the string "a_b_c_α_β" and the label of p_{a} is encoded as the string "a".
    * A word in F-, P-, and R- labels is a multiplication of F-,P-, and R-symbols each raised to a specific power. A word will be represented by a vector w of vectors w[i] = [ l_i, p_i ] where l_i is a vector whose elements represent the labels of an F-,P-, or R-symbol and where p_i is an integer representing the power to which the symbol appears in the word. 
    * For all technical definitions we refer to doi: 10.1090/surv/205.

    The interpretation of the values of the fields of a (pivotal) (braided) fusion category is the following.
    * fusion_ring: a uuid of the Grothendieck ring of the fusion category
    * uuid: UUID1 string that uniquely represents the fusion category. It is independent of any property of the fusion ring and will therefore not change if a property is found to be incorrect, e.g., due to an incorrect property.
    * anyonwiki_code: list of 7 integers r, m, nnsd, i, f, b, p that uniquely identify a fusion category. Here r, m, nnsd, i are the 4 labels of the fusion ring of the category and the f, b, p labels distinguish, respectively, between non-equivalent associators, braided structures, and pivotal structures. Here two braided (resp. pivotal) structures are considered non-equivalent if the braided (resp. pivotal) fusion categories with those structures are not braided (resp. pivotal)-equivalent. The value of b is 0 if the category admits no braiding.
    * f_symbols: JSON dictionary that maps sofs of order 10, representing F-symbols, to either qqb_ids or af_ids representing the exact value of the F-symbol in some specific gauge.
    * r_symbols: JSON dictionary that maps sofs of order 5, representing R-symbols, to qqb_ids or af_ids representing the exact value of the R-symbol in some specific gauge. If the category is not braided, this is an empty JSON object.
    * p_symbols: JSON dictionary that maps sofs of order 1, representing P-symbols, to qqb_ids or af_ids representing the exact value of the P-symbol in some specific gauge. If the category is not pivotal, this is an empty JSON object.
    * quantum_dimensions: vector of qqb_ids [d1,...,dn] where di is the quantum dimension of the i'th simple object  
    * s_matrix: matrix of qqb_ids that represents the S-matrix of the fusion category if the category is a ribbon category
    * twists: vector [θ1,...,θn] of qqb_ids that represent the fractions θi appearing such that the T-matrix obeys T[i,j] = δ_{i,j}exp(2πI*θi)
    * numeric_f_symbols: JSON dictionary that maps sofs of order 10, representing F-symbols, to ctfs representing the numeric value of the F-symbol in some specific gauge.
    * numeric_r_symbols: JSON dictionary that maps sofs of order 5, representing R-symbols, to ctfs representing the numeric value of the R-symbol in some specific gauge.
    * numeric_p_symbols: JSON dictionary that maps sofs of order 1, representing P-symbols, to ctfs representing the numeric value of the P-symbol in some specific gauge.
    * numeric_quantum_dimensions: vector of cfts that represents a numeric approximation to quantum_dimensions 
    * numeric_s_matrix: matrix of ctfs that represents a numeric approximation to s_matrix
    * numeric_twists: vector of ctfs that represents a numeric approximation to the twist_factors
    * base_field: the string "QQBar" if the symbols of the fusion category are expressed using qqb_ids. If the symbols are expressed using af_ids then this is a string of the form "QQ[x]/<a0_a1_..._an>" where the integers a0 to an are the coefficients of the defining polynomial a0 + a1*x + ... + an * x^n of the abstract field.   
    * embedding: if the base field is abstract then a qqb_id that represents the embedding of the generator of the abstract field into QQBar. Otherwise this equals "-1_1__1".
    * minimal_fields: a JSON dictionary that maps the strings "F", "FR", "FP", and "FPR" to the generators (given by qqb_ids) of the mimimal fields for the F-symbols, F-and R-symbols, F-and P-symbols, and F-, P-, and R-symbols respectively. The fields are 
    * in_minimal_field: true if the F-symbols, R-symbols and P-symbols are expressed in the minimal field for all symbols together and false otherwise. 
    * in_unitary_gauge: true if the F-matrices, R-matrices and P-matrices (seen as 1 by 1 matrices of pivotal coeffients) are unitary matrices.
    * names: a JSON dictionary mapping naming conventions to lists of strings of names given using that convention. The conventions at the moment are
      * "quantum_group_like": names associated to quantum groups at level k, such as "[psu(2)_5]_2_2_1", "[so(5)_2]_1_0_1"
      * "group_like": names associated to the theory of finite groups, such as "[Z_2]_1_2_1", "[Rep(D_6)]_2_0_1". Names associated to near-group or group theoretical fusion rings do not belong here but in miscelaneous.
      * "physics": names associated to applications in physics, such as "[Fibonacci]_1_1_1", "[Ising]_2_4_1", "[Potts]_1_0_2".
      * "miscellaneous": names not belonging to another of the above categories such as "[TY(Z_4)]_3_0_1" and "MR_6".
    * texnames: a JSON dictionary mapping naming conventions to lists of strings of names typeset in LaTeX given using that convention. The conventions are the same as for the names field.
    * non_trivial_sub_fusion_cats: list of vectors [ els, uuid ] where els are the simple objects of the parent category that form a subcategory isomorphic to cat[uuid].
    * gauge_split_basis: an array [ I, D ] where I and D are arrays of words in F-,P-, and R-symbols, that represents a gauge-split basis for the fusion category (see: https://doi.org/10.48550/arXiv.2601.20012 for a definition)
    * gauge_split_transform: a vector [ sm, n ] where sm is a matrix representing a sparse array, say M, in CSC format and n is an integer. The matrix M is such that if all F-symbols, R-symbols (if there are any), and P-symbols are sorted lexicographically on their labels and joined in a list in that order, say L, then the elements x_i = ∏_{j}(L[j])^{M[j,i]}
      * form the gauge-invariant part I of the gauge_split_basis for i <= m. Of those elements the ones whose values evaluate to 0 appear first.
      * form the gauge-dependent part D of the gauge_split_basis for i > m 
    * is_pivotal: true if the fusion category is pivotal, false if not
    * is_spherical: true if the fusion category is spherical, false if not
    * is_unitary: true if the fusion category is unitary, false if not. Here by unitary we mean that there exists a gauge in which all F-matrices are unitary matrices and moreover the pivotal structure is such that the quantum dimension of each simple object equals its Frobenius-Perron dimension.
    * is_braided: true if the fusion category is braided, false if not
    * is_ribbon: true if the fusion category is ribbon, false if not
    * is_modular: true if the fusion category is modular, false if not
    * software: JSON dictionary mapping names of fields to a list of reference to software that played a significant role in obtain the data in the way it is represented here. Special field names are
      * "all": when all fields of the category point to the same software
      * "all_other_data": when all other data, besides the data having specific references, points to the same software.
    * references: JSON dictionary mapping names of fields to a list of references to the paper that played a significant role in obtaining the data in the way it is represented here. Special field names are the same as for software. Only papers that have lead to the data as currently represented are included and thus no papers that represent theory that was not directly used, or ,e.g. , data in another format that was not used to obtain current data.
"""

#TODO: add
#   * realizations:
#   * automorphisms
#   
# to the infostring 




function import_fusion_categories( filename::String  )
  jsdict = JSON.parsefile(filename);

  catlist = SkeletalFusioncat[]

  k = keys(jsdict)

  indices = "order" ∈ k ? jsdict["order"] : eachindex(jsdict["data"])
  
  for ind in indices
    
    js = jsdict["data"][ind]

    fr = fawc( fusion_ring_id_to_code(js["fusion_ring"]  ) )

    fs = nothing
      
    rs = nothing

    ps = nothing


    push!( catlist, cat )
  end

  return catlist
  
end



function aw_cat_id( r, m, nnsd, i, j, k, l )
  str = string.( ( r, m, nnsd, i, j, k, l ) )
  "anyonwiki_fcrm_fc__" * join( str, "_" ) 
end

function fusion_ring_id_to_code( id::String )
  c = Parse.( Int, split( id[20:end], "_" ) )

  length(c) ≠ 4 && error("fusion ring should have 4 integers as code. Id most likely invalid") 

  c
end



function jstoqqbfsymbols(jso::JSON.Object{String,Any})
  labels = keys(jso)

  dict = Dict{FL,QQBarElem}( strktotup(l) => from_qqb_id(jso[l]) for l in labels )
  
  
end

# Convert string keys to tuples
function strktotup(s::String)
  Tuple(parse.(Int,(split(s,"_"))))
end

# Convert tuple to string key
function tuptostrk(t)
  join( string.(t), "_" )
end


#┏━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━┓
#┃                             Exporting fusion categories                         ┃
#┗━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━┛

#TODO: the following code only works for cats of which we have the anyonwiki
# codes at the moment


function save_symbols(dir::String,id::String,fs::ASFCS)
  svdir = joinpath(dir,id)

  fn(s::String) = joinpath(svdir, s * ".mrdi" )
  sv(nm::String,data) = Oscar.save(nm, data)
  
  dictkeysfn  = fn( "dict_keys" )
  dictvalsfn  = fn( "dict_vals" )
  polyringfn  = fn( "polyring")
  idictkeysfn = fn( "idict_keys" )
  idictvalsfn = fn( "idict_vals" )
  embeddingfn = fn( "embedding" )

  
  sv( dictkeysfn, labels(fs) )
  sv( dictvalsfn, values(fs) )

  sv( polyringfn, polynomial_ring(fs) )
  
  sv( idictkeysfn, keys(invdict(fs) ) )
  sv( idictvalsfn, values(invdict(fs) ) )

  sv( embeddingfn, embedding(fs) )
  
  d = 
    Dict(
      "dict" => Dict( "keys_fn" => dictkeysfn, "values_fn" =>  dictvalsfn ),
      "polyring" => polyringfn,
      "idict" => Dict( "keys_fn" => idictkeysfn, "values_fn" =>  idictvalsfn ),
      "embedding" => embeddingfn,
      "inminfield" => missing_to_nothing(inminfield(fs)),
      "inunitarygauge" => missing_to_nothing(inunitarygauge(fs))
    )
  write_json( svdir * "symbols.json" )
end

function mult_tab_to_json(mt::Array{Int64,3})
  r = size(mt,1)
  return [ [ [ mt[i, j, k] for k in 1:r ] for j in 1:r ] for i in 1:r ]
end

function missing_to_nothing(x)
  ismissing(x) ? nothing : x
end

function fusion_ring_to_json(fc::SkeletalFusionCat)
  "anyonwiki_"* "fcrm_fr__"*join(string.(anyonwiki_code(fusion_ring(fc))), "_")
end

function fusion_cat_to_json(fc::SkeletalFusionCat)
  
end

#┏━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━┓
#┃                 Exporting and importing F-, R-, and P-symbols                   ┃
#┗━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━┛

# For symbols defined as QQBar elements we will use a similar scheme as
# for fusion rings: we store the QQBar elements in a giant mrdi database
# that gets expanded in individual elements the first time the package
# is loaded. These elements get loaded on demand.
# The filenames of the symbols are determined by the categories they're coming from

function fsfilename(catid::String)::String
  catid*"fsymbols"
end 

function save_F_symbols( catid::String, fs::FSymbols )
  fn = fsfilename(catid)
  
end




# For symbols equal to elements of (possibly embedded) abstract fields
# we will use the following import/export scheme
#
# * each abstract number field has an ID associated to it: asnf_id
# at the moment we only support AbsSimpleNumFieldElems belonging to
# finite extensions of ℚ. Later we will generalize.
# 
# * each set of symbols has an ID: fsymb_id, rsymb_id, psymb_id
# the id contains info on the ring of the stored fusion category and the
# abstract field

# * each embedding into ℂ has an ID: emb_id
# containing info on the ring of the stored fusion category and the qqb
# element of the embedded generator: emb_id

# fusion categories are stored abstractly together with a list of embeddings

# ID of an abstract number field
export asnf_id

function asn_id(F::AbsSimpleNumField)::String
  return polstring(minimal_polynomial(gen(F) ) ) * "__QQ"
end

export

function 

function store_asn(x::AbsSimpleNumFieldElem)
  path = joinpath(@__DIR__, "data", "split_data", "abstractfields", asn_id(x), ".mrdi")  
end

function asn_id(x::AbsSimpleNumFieldElem)::String
  return  "__" * asn_id(parent(x))
end


function save_abs_field_num()

end

