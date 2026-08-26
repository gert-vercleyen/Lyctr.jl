using TensorCategories, FusionRings, JSON, ProgressMeter, Oscar, UUIDs
# Note: fabian doesn't have the new cats yet



# Importing fabians categories 
#allcodes = anyonwiki_keys(7)

#allcats = [ anyonwiki(k...) for k in allcodes ]

# create a list of embeddings from AbsSimpleNumField to QQBar
function tofl( x::AcbFieldElem )
  return Float64(real(x)) + Float64(imag(x)) * im
end

function qqbembim(c)
  ϕ = c.embedding
  r = base_ring(c)

  r == QQ && return QQemb
  
  a = gen(r)
  fltim = tofl( ϕ(a) )
  
  rts = roots( QQBar, minimal_polynomial(a) )
  fltrts = tofl.( AcbField(64).( rts ) )

  mindiff = minimum([ abs2(x - y) for x in fltrts for y in fltrts if y != x  ])

  ( min, ind ) =  findmin( x -> abs2( x - fltim ), fltrts )

  try
    if min > mindiff/2
      println("Wanted: ",fltim)
      println("Closest: ", rts[ind])
      println("Distance: ", abs2(fltim - fltrts[ind] ) )
      println("Minimal distance between roots: ", mindiff)
      println("Accuracy not high enough to guarantee correctness")
    end
  catch e
    print("check embedding ", c)
  end
    
  hom( r, QQBar, rts[ind] )
end

#allembeddings = [ qqbembim(c) for c in allcats ] 

# convert Fabian's F-symbols to my format

function tofl( x::AcbFieldElem )
  return Float64(real(x)) + Float64(imag(x)) * im
end

function fix_F_symbols( cat, emb, p )
  mt = multiplication_table(cat)
  fr = Lyctr.fusion_ring(mt)
  
  formal_fs = F_symbols(fr, field = base_ring(cat), embedding = emb)
  d = F_symbols(cat)
  
  fix_key( v ) = (v[p[1]],v[p[2]],v[p[3]],v[p[4]],v[p[5]],1,1,v[p[6]],1,1) 
  gert_d = Dict( fix_key(v) => val for (v,val) in d  )

  F_symbols( gert_d, polynomial_ring(formal_fs), invdict(formal_fs), emb, mt )
end

function fix_F_symbols_no_perm( cat, emb )
  mt = multiplication_table(cat)
  fr = Lyctr.fusion_ring(mt)
  
  formal_fs = F_symbols(fr, field = base_ring(cat), embedding = emb)
  d = F_symbols(cat)
  
  fix_key( v ) = (v[1],v[2],v[3],v[4],v[1],1,1,v[2],1,1) 
  gert_d = Dict( fix_key(v) => val for (v,val) in d  )


  F_symbols( gert_d, polynomial_ring(formal_fs), invdict(formal_fs), emb, mt )
end


function check_pentagon_equations(cat,perm)
  r     = rank(cat)
  mt    = multiplication_table(cat)
  fsymb = Dict( k[perm] => v for (k,v) in F_symbols(cat) )
  labels = collect(keys(fsymb))
  K = base_ring(cat)


  fs(v) = !haskey(fsymb, v) ? K(0) : fsymb[v]

  # check polynomials for which pent eqn has nonzero LHS
  nl1 = length(labels)
  @showprogress desc="Checking eqns with nonzero lhs" dt=1 for i in 1:nl1
    f, c, d, e, g, l = labels[i]

    function is_match(v)
      return v[3] == l && v[4] == e && v[5] == f
    end

    matches = filter(is_match, labels)
    nl2 = length(matches)
    for j in 1:nl2
      a, b, _, _, _, k = matches[j]

      fac1     = fs([f, c, d, e, g, l])
      fac2     = fs([a, b, l, e, f, k])
      fac3(hh) = fs([a, b, c, g, f, hh])
      fac4(hh) = fs([a, hh, d, e, g, k])
      fac5(hh) = fs([b, c, d, k, hh, l])

      lhs = fac1 * fac2
      rhs = sum( fac3(h) * fac4(h) * fac5(h) for h ∈ 1:r; init = K(0) )
      pol = lhs - rhs

      if pol != K(0)
        println("NONZERO LHS POL")
        @info "a, b, c, d, e, f, g, k, l: " ( a, b, c, d, e, f ,g, k, l )
        @info "lhs" lhs
        @info "rhs" rhs
        return (false,pol)
      end
    end
  end

  # add polynomials for which pent eqn has zero LHS and nonzero RHS
  # This is done by constructing the symmetric tree with non-existent 
  # bottom fusion channel N[f,l,e] and matching the other labels

  fr = fusion_ring(mt)
  sc = nonzero_structure_constants(fr)
  zsc = Tuple{Int64, Int64, Int64}[]

  for i in 1:r, j in 1:r, k in 1:r
    v = (i, j, k)
    if v ∉ sc
      push!(zsc, v)
    end
  end

  @showprogress desc="Checking eqns with zero lhs" dt=1 for n1 in zsc
    f, l, e = n1

    function is_match2(v)::Bool
      return v[3] == f
    end

    function is_match3(v)::Bool
      return v[3] == l
    end

    matches2 = filter(is_match2, sc)
    matches3 = filter(is_match3, sc)
    nl2 = length(matches2)
    nl3 = length(matches3)

    for n2 in matches2, n3 in matches3, k in 1:r, g in 1:r
      a, b, = n2
      c, d, = n3
      fac3(hh) = fs([a, b, c, g, f, hh])
      fac4(hh) = fs([a, hh, d, e, g, k])
      fac5(hh) = fs([b, c, d, k, hh, l])

      pol = sum( fac3(h) * fac4(h) * fac5(h) for h ∈ 1:r; init = K(0) )

      if pol != K(0)
        println("ZERO LHS POL")
        @info "a, b, c, d, e, f, g, k, l: " ( a, b, c, d, e, f ,g, k, l )
        @info "rhs" pol
        return (false,pol)
      end
    end
  end

  return ( true, nothing )
end


# Checking pentagon equations

function group_by_property(f, lis)
  isempty(lis) && return Vector{eltype(lis)}[]

  # Use vector of tuples to preserve order: 
  # (key, ind_of_first_occurrence)
  keys_in_order = Tuple{Any, Int}[]
  groups = Dict{Any, Vector{eltype(lis)}}()

  for (index, el) in enumerate(lis)
    key = f(el)
    # if key is new 
    if !haskey(groups, key)
      # create collection for key 
      groups[key] = eltype(lis)[]
      # store index of key to 
      push!(keys_in_order, (key, index))
    end
    # push element to group whose values evaluate to f(el)
    push!(groups[key], el)
  end

  # Sort keys_in_order by index
  sort!(keys_in_order; by = x -> x[2])

  return [groups[key] for (key, _) in keys_in_order]
end

allcodes = anyonwiki_keys(7)

#allcats = [ anyonwiki(k...) for k in allcodes ]
#nonequivkeys = first.(group_by_property(k -> k[1:5], anyonwiki_keys(7) ) )

#nonequivcats = []
#@showprogress desc="Loading cats" for k in nonequivkeys
#  push!( nonequivcats, anyonwiki(k...) )
#end 

function check_pentagon_equations2(cat,perm;inv = false)
  fsymb =
    if inv
      Dict( k[perm] => v for (k,v) in inverse_F_symbols(cat) )
    else
      Dict( k[perm] => v for (k,v) in F_symbols(cat) )
    end
  
  r = rank(cat)
  mt = multiplication_table(cat)
  n_checked = 0
  K = base_ring(cat)

  getF(v) = haskey(fsymb,v) ? fsymb[v] : throw("KeyNotFound")

  @inbounds for a in 1:r, b in 1:r, c in 1:r, d in 1:r
    for f in 1:r
      mt[a, b, f] == 0 && continue      # a×b -> f
      for l in 1:r
        mt[c, d, l] == 0 && continue     # c×d -> l
        for g in 1:r
          mt[f, c, g] == 0 && continue   # f×c -> g  (g = total of a,b,c)
          for k in 1:r
            mt[b, l, k] == 0 && continue   # b×l -> k (k = total of b,c,d)
            for e in 1:r
              mt[g, d, e] * mt[f, l, e] * mt[a, k, e] == 0 && continue
              # note: we don't check equations where LHS = 0 yet.

              lhs =
                try
                  getF([f, c, d, e, g, l]) * getF([a, b, l, e, f, k])
                catch error
                  return false,  ( a, b, c, d, e, f, g, k, l, "lhs", "KeyNotFound" )
                end

              rhs = K(0)
                try
                  for h in 1:r
                    mt[b, c, h] == 0 && continue   # b×c -> h
                    mt[a, h, g] == 0 && continue   # a×h -> g
                    mt[h, d, k] == 0 && continue   # h×d -> k
                    rhs += getF([a, b, c, g, f, h]) *
                        getF([a, h, d, e, g, k]) *
                        getF([b, c, d, k, h, l])
                  end
                catch error
                    return false,  ( a, b, c, d, e, f, g, k, l, "RHS", "KeyNotFound" )
                end

              if lhs - rhs != 0
                return false, ( a, b, c, d, e, f, g, k, l, lhs, rhs )
              end

            end
          end
        end
      end
    end
  end

  return true, (0,0,0,0,0,0,0,0,0,0,0,0)
end


function group_by_property(f, lis)
  isempty(lis) && return Vector{eltype(lis)}[]

  # Use vector of tuples to preserve order: 
  # (key, ind_of_first_occurrence)
  keys_in_order = Tuple{Any, Int}[]
  groups = Dict{Any, Vector{eltype(lis)}}()

  for (index, el) in enumerate(lis)
    key = f(el)
    # if key is new 
    if !haskey(groups, key)
      # create collection for key 
      groups[key] = eltype(lis)[]
      # store index of key to 
      push!(keys_in_order, (key, index))
    end
    # push element to group whose values evaluate to f(el)
    push!(groups[key], el)
  end

  # Sort keys_in_order by index
  sort!(keys_in_order; by = x -> x[2])

  return [groups[key] for (key, _) in keys_in_order]
end


function formal_F_matrices(flabels)
  grouped_labels = group_by_property(l -> l[1:4], flabels)

  function formal_fmat(lbls)
    # associativity of fusion ring assures size of lbls is square
    m = Int64(sqrt(size(lbls, 1)))
    return reshape(lbls, m, m)
  end

  return formal_fmat.(grouped_labels)
end

function value(fsymb, fmat)
  n = size(fmat,1)
  [ fsymb[fmat[i,j]] for i in 1:n, j in 1:n ]
end


function inverse_F_symbols(cat)
  fs = F_symbols(cat)
  # matrices of labels of F-symbols
  lmats = formal_F_matrices(collect(keys((fs))))
  # matrices of values of F-symbols
  R  = base_ring(cat)
  FF = fraction_field(R)

  Fmats = map( m -> matrix( FF, value(fs, m) ), lmats )

  # since the labels of an F-matrix denote the anyon
  # types the label sets need to be swapped for inverse
  # F symbols to make sense 
  swap_labels(tt) = tt[[1, 2, 3, 4, 6, 5]]
  function swap_labels(m::Matrix)
      [ swap_labels( m[i,j] ) for i in 1:nrows(m), j in 1:ncols(m) ]
  end

  invlmats = swap_labels.(lmats)

  invFmats = inv.(Fmats)

  d = Dict{Vector{Int64},FieldElem}()
  for ind ∈ eachindex(invlmats)
    for i ∈ 1:nrows(invlmats[ind]), j ∈ 1:ncols(invlmats[ind])
      push!(d, invlmats[ind][i,j] => invFmats[ind][i,j] )
    end
  end

  d
end




S_6 = symmetric_group(6)

@showprogress for el in S_6
  println("=== Checking pentagon eqns for label permutation ", el)
   for cat in nonequivcats[2:end] #don't need to check trivial cat 
    println("> ",cat.name)
    check = check_pentagon_equations2(cat,Vector(el),inv=true)
    !(typeof(check)<:Tuple) && throw(cat.name)
    if first(check) == false
      println(cat.name, " : false eqn: ", check[2])
      println("")
      break 
    end 
  end
end




#=
testcases = [ 11, 894 ]

testcats = allcats[testcases]
testembs = allembeddings[testcases]

allfsymbols = [ fix_F_symbols(allcats[i],allembeddings[i]) for i in 1:length(allcats) ]
S1 = symmetric_group(6)

perms = []

maxval = 2

for σ1 ∈ S
  p1 = Vector(σ1,6)
  for i in 1:maxval
    local fs = [ fix_F_symbols( testcats[i], testembs[i], p1 ) for i in 1:maxval ]
    global val = first( check_pentagon_equations(fs[i]) ) 
    val === false && break
  end 
  if val == true
    println((p1), ": TRUE :)")
    push!(perms,(p1))
  else
    println((p1), ": false")
  end
end

perms

=#






# ============= IMPORTING AND CONVERTING SYMBOLS FROM MATHEMATICA ======

# To test the code we only import cats up to maxindex
#=
maxindex = 989

symbols_dir = "/home/gert/Projects/Lyctr.jl/src/data/mathematica/Symbols/"
fns = readdir(symbols_dir)[1:maxindex]
fn_to_code(fn::String) = parse.(Int64,split(fn[9:end-3],"_"))
cat_codes = fn_to_code.(fns)

symboldata = try
    if length(symboldata) < maxindex
      @showprogress desc="loading numeric data" [ [ fn_to_code(fn), include(symbols_dir*fn) ] for fn in fns ]
    else
      symboldata
    end
  catch e
    @showprogress desc="loading numeric data" [ [fn_to_code(fn), include(symbols_dir*fn) ] for fn in fns ]
  end

symboldict = Dict( t[1] => t[2] for t in symboldata )
symbolkeys = collect(keys(symboldict))

function extractfsymbols( data ) 
  fr = fusion_ring_from_data(data)
  formal_fs = F_symbols(fr,field = QQBar, embedding = nothing)
  dict = data["fsymbols"]

  newdict = Dict( fixkey(k) => toqqb(v) for (k,v) in dict )

  F_symbols( newdict, polynomial_ring(formal_fs), invdict(formal_fs), nothing, multiplication_table(fr) )
end

function extractrsymbols( data ) 
  dict = data["rsymbols"]
  isempty(dict) && return nothing
  fr = fusion_ring_from_data(data)
  formal_Rs = R_symbols(fr,field = QQBar, embedding = nothing)

  newdict = Dict( fixkey(k) => toqqb(v) for (k,v) in dict )

  R_symbols( newdict, polynomial_ring(formal_Rs), invdict(formal_Rs), nothing, multiplication_table(fr) )
end

function extractpsymbols( data ) 
  fr = fusion_ring_from_data(data)
  formal_Ps = P_symbols(fr,field = QQBar, embedding = nothing)
  dict = data["psymbols"]

  newdict = Dict( fixkey(k) => toqqb(v) for (k,v) in dict )

  P_symbols( newdict, polynomial_ring(formal_Ps), invdict(formal_Ps), nothing, multiplication_table(fr) )
end

function fusion_ring_from_data( data )
  f_labels = collect(keys(data["fsymbols"]))
  rightvacuumfusion(lab) = lab[3] == 1 
  vacfus = filter( rightvacuumfusion , f_labels )
  structconst = [ ( lab[1], lab[2], lab[4] ) for lab in vacfus  ]

  rank = maximum( l[1] for l in structconst  )

  arr = zeros(Int64,rank,rank,rank)
  for l in structconst
    arr[l[1],l[2],l[3]] = 1
  end
  Lyctr.fusion_ring(arr)
end

function toqqb(tuple)
  R,x = QQ["x"];
  fltval = tuple[2] 

  minpol = R(tuple[1])
  rts = roots( QQBar, minpol )
  l = length(rts)
  fltrts = tofl.( AcbField(64).( rts ) )
  
  ( min, ind ) =  findmin( x -> abs2( x - fltval ), fltrts )

  mindiff = minimum([ abs2(fltrts[i] - fltrts[j]) for i in 1:l for j in (i+1):l ]; init = 1)
  
  if min > mindiff/2
    println("")
    println("Wanted: ",fltval)
    println("Closest: ", rts[ind])
    println("Distance: ", abs2(fltval - fltrts[ind] ) )
    println("Minimal distance between roots: ", mindiff)
    println(tuple)
    println("Accuracy not high enough to guarantee correctness")
    println("-------")
  end

  rts[ind]
end

function tofl( x::AcbFieldElem )
  return Float64(real(x)) + Float64(imag(x)) * im
end

function fixkey(v::Vector{Int64})
  if length(v) == 6
    (v[1],v[2], v[3], v[4], v[5], 1, 1 , v[6], 1, 1 )
  elseif length(v) == 3
    (v[1],v[2],v[3],1,1)
  elseif length(v) == 1
    (v[1],)
  else
    error("Trying to fix labels with wrong length")
  end
end

allfsymbols = 
  try
    if length(allfsymbols) < maxindex
      @showprogress desc="extracting qqb fs" [ code => extractfsymbols(symboldict[code]) for code in symbolkeys ]
    else
      allfsymbols
    end
  catch e
    @showprogress desc="extracting qqb fs" [ code => extractfsymbols(symboldict[code]) for code in symbolkeys ]
  end
fsymboldict = Dict( allfsymbols... )


allrsymbols = 
  try
    if length(allrsymbols) < maxindex
      @showprogress desc="extracting qqb rs" [ code => extractrsymbols(symboldict[code]) for code in symbolkeys ]
    else
      allfsymbols
    end
  catch e
    @showprogress  desc="extracting qqb rs" [ code => extractrsymbols(symboldict[code]) for code in symbolkeys ]
  end
rsymboldict = Dict( allrsymbols... )


allpsymbols = 
  try
    if length(allpsymbols) < maxindex
      @showprogress desc="extracting qqb ps" [ code => extractpsymbols(symboldict[code]) for code in symbolkeys ]
    else
      allfsymbols
    end
  catch e
    @showprogress desc="extracting qqb ps" [ code => extractpsymbols(symboldict[code]) for code in symbolkeys ]
  end

psymboldict = Dict( allpsymbols... )
=#
function write_json(filename::String, data::Dict)
  open(filename, "w") do f
    return JSON.json(f, data; pretty = true, inline_limit = 10)
  end
end

#=
inviddict =
  if "id_dict.json" ∈ readdir("/home/gert/Projects/Lyctr.jl/src/data/fusion_categories/")
    println("Using saved uuids")
    JSON.parsefile("/home/gert/Projects/Lyctr.jl/src/data/fusion_categories/id_dict.json")
  else
    println("WARNING: Creating new uuidtocodedict!!!")
    uuidtocodedict = Dict(string(UUIDs.uuid1())=>c for c in cat_codes)
    write_json("/home/gert/Projects/Lyctr.jl/src/data/fusion_categories/id_dict.json",uuidtocodedict)
    uuidtocodedict
  end

iddict = Dict( v => k for (k,v) in inviddict )

catid(c) = iddict[c]


qqb_id = Lyctr.qqb_id

qqbd = Lyctr.qqb_dict

datadir = "/home/gert/Projects/Lyctr.jl/src/data/"
println("importing ids")
newids  = Oscar.load(joinpath( datadir, "qqb_ids.mrdi"))
println("importing vals")
newvals = Oscar.load(joinpath( datadir, "qqb_vals.mrdi"))
println("Finished")


@showprogress desc="adding new vals" for i in 1:length(newids)
  qqbd[newids[i]] = newvals[i]
end


invqqbd = Dict( v => k for (k,v) in qqbd )
=#
#=
function qqbid(v)
  if haskey(invqqbd,v)
    invqqbd[v]
  else
    id = qqb_id(v)
    invqqbd[v] = id
    qqbd[id] = v
    Oscar.save(joinpath(datadir, "split_number_data", invqqbd[v]*".mrdi" ), v )
    id
  end
end

function exportfsymbols()
  allfsymbolsd = Dict{String,Dict{String,String}}()
  @showprogress desc="exporting fsymbols" for code in symbolkeys
    fs = fsymboldict[code]
    fd = dict(fs)
    
    # populate qqb_dict
    exportedfdict = Dict{String,String}()
    for (k,v) in fd
      key  = join( string.(k), "_" )
      qqid = qqbid(v)
      
      exportedfdict[key] = qqid  
    end
    allfsymbolsd[ catid(code) ] = exportedfdict
  end

  data = Dict(
    "data" => allfsymbolsd,
    "order" => [ catid(code) for code in symbolkeys ],
    "info" =>
    """
      This file contains data on F-symbols for multiplicity-free fusion categories.
      The
        * \"data\" field contains a dictionary mapping identifiers of fusion categories (fc_ids) to dictionaries (fc_dicts).
          Here 
          * fc_ids are the identifiers for the categories consist of several letters denoting the type of classification, separated via a __ by an identifier for the specific  category. The current lettercodes for classifications are { fcrmfc, which stands for [f]ull [c]lassification of [r]ank and [m]ultiplicity [f]usion [c]ategory which is also known as the anyonwiki code, and is followed by 7 integers }. 
          * fc_dicts are dictionaries mapping (F_labels) to (qqb_ids)
            Here
            * F_labels are strings encoding the 10 indices that label an F-symbols
            * qqb_ids are strings that uniquely encode algebraic numbers.
    """
  )
  
  #export F-symbols
  write_json(
    joinpath( datadir, "fusion_categories/f_symbols_qqb.json" ),
    data
  )
  
  #export qqb_dict
  #Oscar.save( joinpath( datadir, "qqb_ids.mrdi"), collect(keys(qqbd)) )
  #Oscar.save( joinpath( datadir, "qqb_vals.mrdi"), collect(values(qqbd)) )
end


function exportrsymbols()
  allrsymbolsd = Dict{String,Dict{String,String}}()
  @showprogress desc = "exporting rsymbols" for code in symbolkeys
    rs = rsymboldict[code]

    isnothing(rs) && continue

    rd = dict(rs)
    # populate qqb_dict
    exportedrdict = Dict{String,String}()
    for (k,v) in rd
      key  = join( string.(k), "_" )
      qqid = qqbid(v)
      exportedrdict[key] = qqid
    end
    allrsymbolsd[ catid(code) ] = exportedrdict
  end


  data = Dict(
    "data" => allrsymbolsd,
    "order" => [ catid(code) for code in symbolkeys ], 
    "info" =>
    """
      This file contains data on R-symbols for multiplicity-free fusion categories.
      The
        * \"data\" field contains a dictionary mapping identifiers of fusion categories (fc_ids) to dictionaries (fc_dicts).
          Here 
          * fc_ids are the identifiers for the categories consist of several letters denoting the type of classification, separated via a __ by an identifier for the specific  category. The current lettercodes for classifications are { fcrmfc, which stands for [f]ull [c]lassification of [r]ank and [m]ultiplicity [f]usion [c]ategory which is also known as the anyonwiki code, and is followed by 7 integers }. 
          * fc_dicts are  dictionaries mapping (R_labels) to (qqb_ids). If the category is not braided, its corresponding dictionary is empty. If no braiding information is available, no dict is stored, i.e. the fc_id points to nothing.
            Here
            * R_labels are strings encoding the 5 indices that label an R-symbol
            * qqb_ids are strings that uniquely encode algebraic numbers.
    """
  )
  
  #export R-symbols
  write_json(
    joinpath( datadir, "fusion_categories/r_symbols_qqb.json" ),
    data
  )
  
  #export qqb_dict
end

function exportpsymbols()
  allpsymbolsd = Dict{String,Dict{String,String}}()
  @showprogress desc = "exporting psymbols" for code in symbolkeys
    ps = psymboldict[code]
    pd = dict(ps)
    # populate qqb_dict
    exportedpdict = Dict{String,String}()
    for (k,v) in pd
      key  = join( string.(k), "_" )
      qqid = qqbid(v)
      exportedpdict[key] = qqid
        
    end
    allpsymbolsd[ catid(code) ] = exportedpdict
  end
  data = Dict(
    "data" => allpsymbolsd,
    "order" => [ catid(code) for code in symbolkeys ],
    "info" =>
    """
      This file contains data on P-symbols (or pivotal coefficients) for multiplicity-free fusion categories.
      The
        * \"data\" field contains a dictionary mapping identifiers of fusion categories (fc_ids) to dictionaries (fc_dicts).
          Here 
          * fc_ids are the identifiers for the categories consist of several letters denoting the type of classification, separated via a __ by an identifier for the specific  category. The current lettercodes for classifications are { fcrmfc, which stands for [f]ull [c]lassification of [r]ank and [m]ultiplicity [f]usion [c]ategory which is also known as the anyonwiki code, and is followed by 7 integers }. 
          * fc_dicts are dictionaries mapping (P_labels) to (qqb_ids)
            Here
            * F_labels are strings encoding the index that label a P-symbol
            * qqb_ids are strings that uniquely encode algebraic numbers.
    """
  )
  
  #export P-symbols
  write_json(
    joinpath( datadir, "fusion_categories/p_symbols_qqb.json" ),
    data
  )
end

#exportfsymbols()
#exportrsymbols()
#exportpsymbols()

#=
println("Loading data")
inv_id_dict = JSON.parsefile("/home/gert/Projects/Lyctr.jl/src/data/fusion_categories/id_dict.json")
tocode(s::String) = parse.(Int64,split(s,"_")) 
id_dict = Dict( v => k for (k,v) in id_dict  )

all_fsymbols = JSON.parsefile("/home/gert/Projects/Lyctr.jl/src/data/fusion_categories/f_symbols_qqb.json",Dict{String,Any})["data"]
#all_psymbols = JSON.parsefile("/home/gert/Projects/Lyctr.jl/src/data/fusion_categories/p_symbols_qqb.json",Dict{String,Any})["data"]
#all_rsymbols = JSON.parsefile("/home/gert/Projects/Lyctr.jl/src/data/fusion_categories/r_symbols_qqb.json",Dict{String,Any})["data"]
#properties = JSON.parsefile("/home/gert/Projects/Lyctr.jl/src/data/fusion_categories/cat_props.json",Dict{String,Any})
all_qdims  = JSON.parsefile("/home/gert/Projects/Lyctr.jl/src/data/fusion_categories/quantum_dims.json",Dict{String,Any})
all_s_mats = JSON.parsefile("/home/gert/Projects/Lyctr.jl/src/data/fusion_categories/s_matrices.json",Dict{String,Any})
all_twists = JSON.parsefile("/home/gert/Projects/Lyctr.jl/src/data/fusion_categories/spins.json",Dict{String,Any})

println("Loaded all data")

_nonames = Dict(
  "quantum_group_like" => missing,
  "group_like"         => missing,
  "physics"            => missing,
  "miscellaneous"      => missing,
)

function mscnames(v::Vector)
  return Dict(
    "quantum_group_like" => missing,
    "group_like"         => missing,
    "physics"            => missing,
    "miscellaneous"      => string.(v),
  )
end
=#

function fus_cat_to_dict(id::String)
  cat_code = inv_id_dict[id]
  str_code = join(string.(cat_code),"_")
  fr_code  = cat_code[1:4]
  fr       = fawc(fr_code...)
  fr_id    = uuid(fr)

  props = properties[str_code]
  is_braided   = haskey( all_rsymbols, id )
  is_unitary   = props["unitary"]
  is_spherical = props["spherical"]
  is_ribbon    = props["ribbon"]
  is_modular   = props["modular"]
  
  fdict = all_fsymbols[id]
  pdict = all_psymbols[id]
  rdict = is_braided ? all_rsymbols[id] : Dict{String,String}()

  
  function from_qqb_id(s::String)
    if haskey(qqbd, s)
      return qqbd[s]
    else
      fn = joinpath(datadir, "split_number_data", s*".mrdi")

      if isfile(fn)
        val = Oscar.load(fn)
        qqbd[s] = val
        return val
      else
        spl = split(s, "_")
        l   = length(spl)
        n   = parse(Int64, last(spl))
        cfs = ZZ.(parse.(Int64, spl[1:(end - 2)]))

        qqbval = sort(roots(QQBar, Oscar.polynomial(ZZ, cfs)); by = Lyctr.root_sort_crit)[n]

        path = joinpath(datadir, "split_number_data", s*".mrdi")
        Oscar.save(path, qqbval)

        return qqbval
      end
    end
  end

  function qqb_id_to_flt(s::String)
    x = AcbField(64)(from_qqb_id(s))
    return [ Float64(real(x)), Float64(imag(x)) ]
  end
    
  function to_fl(dict)
    Dict{String,Vector{Float64}}( k => qqb_id_to_flt(v) for (k,v) in dict )
  end

  nfdict = to_fl(fdict)
  npdict = to_fl(pdict)
  nrdict = is_braided ? all_rsymbols[id] : Dict{String,Vector{Float64}}()

  qdims = all_qdims[str_code]
  
  nqdims = qqb_id_to_flt.(qdims)

  smat  =
    if is_modular
      [ all_s_mats[str_code][i][j] for i in 1:rank(fr), j in 1:rank(fr) ] 
    else
      nothing
    end
  
  nsmat =
    if is_modular
      [ qqb_id_to_flt(smat[i,j]) for i in 1:rank(fr), j in 1:rank(fr) ]
    else
      nothing
    end
  
  twists = 
    if is_modular
      all_twists[str_code] 
    else
      nothing
    end

  ntwists = 
    if is_modular
      qqb_id_to_flt.(twists)
    else
      nothing
    end
  
  fr_texnames  = FusionRings.tex_names(fr)
  function to_cat_tex_name(s::String)
    i, j, k = cat_code[5:7]
    "["*s*"]_{$i,$j,$k}" 
  end
  texnames = mscnames(to_cat_tex_name.(fr_texnames["miscellaneous"]))
    
  Dict(
    "fusion_ring"                   => fr_id,
    "uuid"                          => id, 
    "anyonwiki_code"                => cat_code,
    "f_symbols"                     => fdict,
    "r_symbols"                     => rdict,
    "p_symbols"                     => pdict,
    "quantum_dimensions"            => qdims,
    "s_matrix"                      => smat,
    "twists"                        => twists,
    "numeric_f_symbols"             => nfdict,
    "numeric_r_symbols"             => nrdict,
    "numeric_p_symbols"             => npdict,
    "numeric_quantum_dimensions"    => nqdims,
    "numeric_s_matrix"              => nsmat,
    "numeric_twists"                => ntwists,
    "base_field"                    => "QQBar",
    "embedding"                     => nothing,
    "minimal_fields"                => nothing,
    "in_minimal_field"              => nothing,
    "in_unitary_gauge"              => is_unitary,
    "names"                         => _nonames,
    "tex_names"                     => texnames,
    "non_trivial_sub_fusion_cats"   => nothing,
    "gauge_split_transform"         => nothing,
    "gauge_split_basis"             => nothing,
    "realizations"                  => Dict{String,Any}(),
    "is_pivotal"                    => true,
    "is_spherical"                  => is_spherical,
    "is_unitary"                    => is_unitary,
    "is_braided"                    => is_braided,
    "is_modular"                    => is_modular,
    "references"                    => Dict("All" => "arXiv:2405.20075"),
    "software"                      => Dict("All" => "https://doi.org/10.5281/zenodo.10686859")
  )
end

ids = collect(values(id_dict))


@showprogress for i in 1:989
  fn = datadir*"fusion_categories/allcatdata/"*"cat_$i.json"
  FusionRings.write_json( fn, fus_cat_to_dict(ids[i]) )
end

dicts = [ id => fus_cat_to_dict(id) for id in ids ]

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

allcats = Dict(
  "order" => ids,
  "info" => "List of multiplicity-free complex fusion categories.\n\n" * fusion_cat_info_text,
  "data" => Dict(dicts...)
)

write_json( datadir * "MultFreeFusionCategories.json" , allcats, )

#=
datadir =  "/home/gert/Projects/Lyctr.jl/src/data/fusion_categories/"
fdata = JSON.parsefile(joinpath(datadir,"unitary_f_symbols_qqb.json"))
rdata = JSON.parsefile(joinpath(datadir,"unitary_r_symbols_qqb.json"))
pdata = JSON.parsefile(joinpath(datadir,"unitary_p_symbols_qqb.json"))

qqbdict =
  begin
    ks = Oscar.load("src/data/qqb_ids.mrdi")
    vls = Oscar.load("src/data/qqb_vals.mrdi")
    Dict( (ks[i], vls[i]) for i in 1:length(ks) )
  end

# Convert string keys to tuples
function strktotup(s::String)
  Tuple(parse.(Int,(split(s,"_"))))
end

function tuptostrk(t)
  join( string.(t), "_" )
end

function fusion_ring_from_fsymbols( fdict )
  c = collect(keys(fdict))
  f_labels = strktotup.(c)
  rightvacuumfusion(lab) = lab[3] == 1 
  vacfus = filter( rightvacuumfusion , f_labels )
  structconst = [ ( lab[1], lab[2], lab[4] ) for lab in vacfus  ]

  rank = maximum( l[1] for l in structconst  )


  arr = zeros(Int64,rank,rank,rank)
  for l in structconst
    arr[l[1],l[2],l[3]] = 1
  end
  Lyctr.replace_by_known(Lyctr.fusion_ring(arr))
end


function fus_cat_from_data(catid)
  fdict = fdata["data"][catid]

  rdict = try
    rdata["data"][catid]
  catch e
    nothing
  end
  
  pdict = pdata["data"][catid]

  # FUSION RING
  fr = fusion_ring_from_fsymbols(fdict)
  mt = multiplication_table(fr)

  function datatodict( sdict )
    Dict( strktotup(k) => qqbdict[sdict[k]] for k in keys(sdict) )
  end

  # ID 
  id = catid

  # F-, R-, P-, and G-SYMBOLS

  # construct F-symbols
  fd      = datatodict(fdict)
  nf      = length(fd)
  flabels = sort(collect(keys(fd)))
  Rf, f   = polynomial_ring(QQBar, :ℱ => 1:nf )
  invdict = Dict( f[i] => flabels[i] for i in 1:nf )
  
  fsymb = F_symbols(
    datatodict(fdict),
    Rf,
    invdict,
    nothing,
    mt,
    inunitarygauge=true
  )

  # construct R-symbols
  rsymb = if isnothing(rdict)
    nothing
    else 
      rd      = datatodict(rdict)
      nr      = length(rd)
      rlabels = sort(collect(keys(rd)))
      Rr, r   = polynomial_ring(QQBar, :ℛ => 1:nr )
      invdict = Dict( r[i] => rlabels[i] for i in 1:nr )

      R_symbols(
      datatodict(rdict),
      Rr,
      invdict,
      nothing,
      mt,
      inunitarygauge=true
    )
  end
  

  # construct P-symbols
  pd      = datatodict(pdict)
  np      = length(pd)
  plabels = sort(collect(keys(pd)))
  Rp, p   = polynomial_ring(QQBar, :𝒫 => 1:np )
  invdict = Dict( p[i] => plabels[i] for i in 1:np )

  psymb = P_symbols(
    datatodict(pdict),
    Rp,
    invdict,
    nothing,
    mt,
    inunitarygauge=true
  )

  # construct G-symbols
  gsymb = G_symbols(fr, field = QQBar, embedding = nothing)

  # FIELD
  K = QQBar

  # EMBEDDING
  emb = nothing

  # INMINFIELD
  inminfield = missing
  
  # MINIMALFIELDS
  minfields = missing

  skeletal_fusion_cat(
    fr, id, fsymb, K, emb,
    rsymbols = rsymb, psymbols = psymb, inunitarygauge = true
  )
end

#cats = @showprogress [ fus_cat_from_data(id) for id in fdata["order"] ]

qqbinvdict = Dict( v => k for (k,v) in qqbdict )

function fix_fusion_cat_id( id::String )
  if id[1:5] == "crmfc"
    "anyonwiki_fcrm_fc__" * id[6:end]
  else
    id
  end
end

function fusion_ring_id( code )
  "anyonwiki_fcrm_fr__"* join(string.(code),"_")
end

function fusion_cat_to_dict(cat::SkeletalFusionCat)
  fr = fusion_ring(cat)
  c = Lyctr.anyonwiki_code(fr)

  fdict = dict(F_symbols(cat))
  fsymb = Dict( tuptostrk(k) => qqbinvdict[v] for (k,v) in fdict )

  rs = R_symbols(cat)
  rsymb =
    if ismissing(rs) || isnothing(rs)
      nothing
    else
      Dict( tuptostrk(k) => qqbinvdict[v] for (k,v) in dict(rs) )
    end
  isbr = !isnothing(rs)

  pdict = dict(P_symbols(cat))
  psymb = Dict( tuptostrk(k) => qqbinvdict[v] for (k,v) in pdict )
    
  data =
    Dict(
      "fusion_ring"                   => fusion_ring_id(Lyctr.anyonwiki_code(fr)),
      "id"                            => fix_fusion_cat_id(fusion_cat_id(cat)),
      "f_symbols"                     => fsymb,
      "r_symbols"                     => rsymb,
      "p_symbols"                     => psymb,
      "base_field"                    => "QQBar",
      "embedding"                     => nothing,
      "in_minimal_field"              => nothing,
      "minimal_fields"                => nothing,
      "in_unitary_gauge"              => nothing,
      "names"                         => nothing,
      "tex_names"                     => nothing,
      "non_trivial_sub_fusion_cats"   => nothing,
      "inverse_gauge_split_transform" => nothing,
      "gauge_split_basis"             => nothing,
      "realizations"                  => nothing,
      "is_pivotal"                    => true,
      "is_unitary"                    => nothing,
      "is_braided"                    => isbr,
      "is_modular"                    => nothing,
      "references"                    => Dict("All" => "arXiv:2405.20075"),
      "software"                      => Dict("All" => "https://doi.org/10.5281/zenodo.10686859")
    )
end

fusion_cat_info_text = """
    The following conventions are used in explaining the values of the dictionary:
    * While the data describes a (pivotal) (braided) fusion category, we will always refer to this category as a fusion category to save space.
    * A qqb_id is a string that uniquely describes an algebraic number. It is formatted as a list of n integers a0, ..., an, separated by underscores, followed by a double underscore, followed by an integer i: "a0_a1_..._an__i". Here a0 to an are the coefficients of a polynomial a0 + a1*x + ... + an*x^n and i denotes the i'th root of that polynomial. The indexing of roots of the polynomial takes the real roots first, in increasing order. Then come complex conjugate pairs of roots, sorted first by increasing real part and second by increasing complex part.
    * A af_id is a string that describes an element of an abstract field. It is represented as an array of n integers [a0, ..., an]. If x is the generator of the abstract field (stored under the key base_field as part of the data of a fusion category), then the vector represents the element a0 + a1*x + ... + an*x^n where n is the order of the field minus 1.
    * A ctf, or complex tuple of floats, is a representation of a complex floating point number a + i b by vector [ a, b ] where a and b are real floating point numbers.
    * If A is a JSON array then A[i] is its i'th element.
    * The mother category is the current category being represented by the JSON dictionary. When talking about subcategories, the mother category represents the category of which we're interested in its subcategories.
    * Every fusion category is uniquely identified by a uuid. By cat[cat_uuid] we mean the fusion category with uuid equal to cat_uuid
    * All stored fusion categories have a fixed order of their simple objects which equals the order of the elements of the fusion ring that is referenced. Therefore every simple object of a fusion category can be represented by a positive integer from 1 to r = rank(cat) we will often use elements and indices representing those elements interchangeably. By the i'th element of the cat[uuid] we mean the i'th element for the stored fusion category unless stated otherwise.
    * A value of null for any field means the data is missing unless states otherwise. In particular, it does not imply that the data doesn't exist.
    * A sof, or string of labels, is a string of integers seperated by underscores, i.e. of the form "a1_a2_..._an" where each ai is an integer. Its order is by definition the number of integers in the string.
    * By a P-symbol we mean a pivotal coefficient. (see "G. Vercleyen, On Low-Rank Multiplicity-Free Fusion Categories. National University of Ireland, Maynooth (Ireland), 2024." for the definition) 
    * The labels of the F-symbols, R-symbols, and P-symbols follow the conventions of "G. Vercleyen, On Low-Rank Multiplicity-Free Fusion Categories. National University of Ireland, Maynooth (Ireland), 2024." Where the labels of the symbol [F^{a,b,c}_d]^{(e,α,β)}_{(f,γ,δ)} are encoded as the string "a_b_c_d_e_α_β_f_γ_δ", the labels of [R^{a,b}_c]^{α}_{β} are encoded as the string "a_b_c_α_β" and the label of p_{a} is encoded as the string "a".
    * For all technical definitions we refer to doi: 10.1090/surv/205.

    The interpretation of the values of the fields of a (pivotal) (braided) fusion category is the following.
    * fusion_ring: a uuid of the Grothendieck ring of the fusion category
    * uuid: UUID1 string that uniquely represents the fusion category. It is independent of any property of the fusion ring and will therefore not change if a property is found to be incorrect, e.g., due to an incorrect property.
    * anyonwiki_code: list of 7 integers r, m, nnsd, i, f, b, p that uniquely identify a fusion category. Here r, m, nnsd, i are the 4 labels of the fusion ring of the category and the f, b, p labels distinguish, respectively, between non-equivalent associators, braided structures, and pivotal structures. Here two braided (resp. pivotal) structures are considered non-equivalent if the braided (resp. pivotal) fusion categories with those structures are not braided (resp. pivotal)-equivalent. The value of b is 0 if the category admits no braiding.
     * f_symbols: JSON dictionary that maps sofs of order 10, representing F-symbols, to either qqb_ids or af_ids representing the exact value of the F-symbol in some specific gauge.
    * r_symbols: JSON dictionary that maps sofs of order 5, representing R-symbols, to qqb_ids or af_ids representing the exact value of the R-symbol in some specific gauge. If the category is not braided, this is an empty JSON object.
    * p_symbols: JSON dictionary that maps sofs of order 1, representing P-symbols, to qqb_ids or af_ids representing the exact value of the P-symbol in some specific gauge. If the category is not pivotal, this is an empty JSON object.
    * numeric_f_symbols: JSON dictionary that maps sofs of order 10, representing F-symbols, to ctfs representing the numeric value of the F-symbol in some specific gauge.
    * numeric_r_symbols: JSON dictionary that maps sofs of order 5, representing R-symbols, to ctfs representing the numeric value of the R-symbol in some specific gauge.
    * numeric_p_symbols: JSON dictionary that maps sofs of order 1, representing P-symbols, to ctfs representing the numeric value of the P-symbol in some specific gauge.
    * base_field: the string "QQBar" if the symbols of the fusion category are expressed using qqb_ids. If the symbols are expressed using af_ids then this is a string of the form "QQ[x]/<a0_a1_..._an>" where the integers a0 to an are the coefficients of the defining polynomial a0 + a1*x + ... + an * x^n of the abstract field.   
    * embedding: if the base field is abstract then a qqb_id that represents the embedding of the generator of the abstract field into QQBar. Otherwise this equals "-1_1__1".
    * minimal_fields: a JSON dictionary that maps the strings "F", "FR", "FP", and "FPR" to the generators (given by qqb_ids) of the mimimal fields for the F-symbols, F-and R-symbols, F-and P-symbols, and F-, P-, and R-symbols respectively. The fields are 
    * in_unitary_gauge: true if the F-matrices, R-matrices and P-matrices (seen as 1 by 1 matrices of pivotal coeffients) are unitary matrices.
    * names: a JSON dictionary mapping naming conventions to lists of strings of names given using that convention. The conventions at the moment are
      * "quantum_group_like": names associated to quantum groups at level k, such as "[psu(2)_5]_2_2_1", "[so(5)_2]_1_0_1"
      * "group_like": names associated to the theory of finite groups, such as "[Z_2]_1_2_1", "[Rep(D_6)]_2_0_1". Names associated to near-group or group theoretical fusion rings do not belong here but in miscelaneous.
      * "physics": names associated to applications in physics, such as "[Fibonacci]_1_1_1", "[Ising]_2_4_1", "[Potts]_1_0_2".
      * "miscellaneous": names not belonging to another of the above categories such as "[TY(Z_4)]_3_0_1" and "MR_6".
    * texnames: a JSON dictionary mapping naming conventions to lists of strings of names typeset in LaTeX given using that convention. The conventions are the same as for the names field.
    * gauge_split_transform: 
    * gauge_split_basis:
    * is_pivotal: true if the fusion category is pivotal, false if not
    * is_spherical: true if the fusion category is spherical, false if not
    * is_unitary: true if the fusion category is unitary, false if not. Here by unitary we mean that there exists a gauge in which all F-matrices are unitary matrices and moreover the pivotal structure is such that the quantum dimension of each simple object equals its Frobenius-Perron dimension.
    * is_braided: true if the fusion category is braided, false if not
    * is_ribbon: true if the fusion category is ribbon, false if not
    * is_modular: true if the fusion category is modular, false if not
"""

function write_json(filename::String, data::Dict)
  open(filename, "w") do f
    return JSON.json(f, data; pretty = true, inline_limit = 10)
  end
end

function export_fusion_cat(filename::String, cat::SkeletalFusionCat)
  dict = fusion_cat_to_dict(cat::SkeletalFusionCat)
  dict["info"] =
    "Skeletal (Pivotal) (Braided) Fusion Category Determined by F (R) (and P) symbols.\n" *  fusioncatfieldinfostring

  write_json(filename,dict)
end

function export_fusion_cats(filename::String,cats::Vector{SkeletalFusionCat})

  data =
    Dict( fusion_cat_id(cat) => fusion_cat_to_dict(cat) for cat in cats)

  order = fusion_cat_id.(cats)

  infostring =
    "Skeletal (Pivotal) (Braided) Fusion Categories Determined by F (R) (and P) symbols.\n" *
    "The field \"order\" contains the ids of the fusion categories in the original order they were exported.\n"*
    "The field \"data\" contains a dictionary mapping ids of fusion categories to dictionaries with their data.\n"*
    fusioncatfieldinfostring
  
  d = Dict(
    "data" => data,
    "order" => order,
    "info" => infostring
  )
  write_json(filename, d )
end

fusion_cat_to_dict(cats[234])

export_fusion_cats("/home/gert/Projects/Lyctr.jl/src/data/fusion_categories/fusion_categories_ug1.json", cats)
=#



#=
for catid in fdata["order"]
  
end
=#

#=

# polynomials over QQ
R, x = polynomial_ring(QQ,"𝑥")

# === CASE 1: combination of embedded fields over same base field === 

# example fields and embeddings to experiment with
# I'm assuming the base fields of F & M are QQ for simplicity
mp1    = x^3 + x + 1
F, a   = number_field( mp1, "a" )
f      = hom( F, QQBar, roots(QQBar,mp1)[2] )

mp2    = cyclotomic_polynomial(5)
M, b   = number_field( mp2, "b" )
QQAb, z = abelian_closure(QQ)
m      = hom( M, QQBar, QQBar( z(5) ) )

# number field generated by generators of F and M
FM, fm = number_field( base_field(F), [ f(a), m(b) ] )

# FM might not be a simple extension?
sFM, sfm = absolute_simple_field(FM)

# sFM can be quite ugly so we simplify it
K, simp = simplify( sFM )

# Now we want to have
# 1. the embedding of elements of K into QQBar
# 2. homs that map elements of F, M into K

# 1. the embedding of K into ℂ
k = hom( K, QQBar, (fm ∘ sfm ∘ simp)(gen(K)) )

# 2. the morphisms from F1 & F2 into S should be
h_FK = hom( F, K, (inv(k) ∘ f )(a) )
h_MK = hom( M, K, (inv(k) ∘ m )(b) )

# but this gives errors :(
# note: computing ϕS takes some time. We might be using number field  
# functionality in a bad way...
#


# === CASE 2: combination of embedded fields F1, F2 where
# base_field(Fd) == F1 


mp1    = x^3 + x + 1
F1, a1 = number_field( mp1, "𝑎₁" )
ϕ1     = hom( F1, QQBar, roots(QQBar,mp1)[2] )

mp2    = change_base_ring(F1,cyclotomic_polynomial(6))
F2, a2 = number_field( mp2, "𝑎₂" )
S, ϕS  = absolute_simple_field(F2)

QQab, ζ = abelian_closure(QQ)
#ϕ2     = hom( F2, QQBar, QQBar(ζ(6)) )


=#
=#
