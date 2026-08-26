#┏━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━┓
#┃                              gauge_symmetries                                   ┃
#┗━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━┛

function gauge_symmetries(symb::FSymbols; symbol=:𝒢 )
  fr = fusion_ring(symb)

  l = labels(symb)

  function φ(tjs)
      value( symb, tjs )
  end
  
end




#┏━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━┓
#┃                            gauge-split_transform                                ┃
#┗━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━┛
export gauge_split_transform

function gauge_split_transform( fr::FusionRing; zeros::Vector{FL} = [], includeonly::Vector{String} = ["FSymbols","PSymbols","RSymbols"] )
  multiplicity(fr) > 1 && error("gauge_split_transform is only implemented for multiplicity-free fusion rings")
  includeonly ⊈ ["FSymbols","PSymbols","RSymbols"] && error("option includeonly should be a subset of [\"FSymbols\",\"PSymbols\",\"RSymbols\"]")

  lio = length(includeonly)
  rank(fr) == 1 && return ( identity_matrix(SMat, ZZ, lio ), lio )
  
  fs = "FSymbols" ∈ includeonly ? dict(FSymbols(fr, triv_vac = false)) : [] 
  fl = isempty(fs) ? [] : sort(collect(keys(fs)))
  fv = isempty(fs) ? [] : [ fs[l] for l in fl ]
  
  rs = "RSymbols" ∈ includeonly ? dict(RSymbols(fr, triv_vac = false)) : [] 
  rl = isempty(rs) ? [] : sort(collect(keys(rs)))
  rv = isempty(rs) ? [] : [ rs[l] for l in rl ]

  ps = "PSymbols" ∈ includeonly ? dict(PSymbols(fr, triv_vac = false)) : []
  pl = isempty(ps) ? [] : sort(collect(keys(ps)))
  pv = isempty(ps) ? [] : [ ps[l] for l in pl ]

  labels  = vcat(fl,rl,pl)
  values  = vcat(fv,rv,pv)
  nlabels = length(labels)

  zeropos    = [ x -> findfirst(x, labels) for x in zeros ]
  nonzeropos = filter( i -> i ∉ zeropos, collect(1:nlabels) )
  # permutation to move all zero elements to front
  σ = invperm(vcat(zeropos,nonzeropos))

  # set up the gauge matrix and transpose
  gaugemat = gauge_matrix(fr,zeros,includeonly=includeonly)

  # we want the kernel of the transposed gauge matrix
  transfmat, upper_triang = hnf_with_transform(tgaugemat)
  
  # rank of gaugemat
  r = nnz(upper_triang)
  # num cols of transformation matrix
  nc = number_of_columns(transfmat)
  
  # add the zero positions
  fulltransfmat =
    if isempty(zeropos)
      transfmat
    else
      lzp = length(zeropos)
      ds = d_sum(
        identity_matrix( SMat, ZZ, lzp ),
        transpose(transfmat)
      )
      sm = sparse_matrix(ZZ, 0, lzp)
      for i in 1:lzp
        push!(sm, ds[σ[i]] )
      end
    end

  return (
    fulltransfmat
    ,
    nc - r + length(zeropos)
  )
  
end

function gauge_matrix(fr::FusionRing; zeros::Vector{FL} = [], includeonly::Vector{String} = ["FSymbols","PSymbols","RSymbols"] ) 
  multiplicity(ring) ≠ 1 && error("gauge_matrix is only defined for multiplicity_free fusion rings")
  
  lio = length(includeonly)
  rank(fr) == 1 && return ( identity_matrix(SMat, ZZ, lio ), lio )
  
  fs = "FSymbols" ∈ includeonly ? dict(FSymbols(fr, triv_vac = false)) : [] 
  fl = isempty(fs) ? [] : sort(collect(keys(fs)))
  fv = isempty(fs) ? [] : [ fs[l] for l in fl ]
  
  rs = "RSymbols" ∈ includeonly ? dict(RSymbols(fr, triv_vac = false)) : [] 
  rl = isempty(rs) ? [] : sort(collect(keys(rs)))
  rv = isempty(rs) ? [] : [ rs[l] for l in rl ]

  ps = "PSymbols" ∈ includeonly ? dict(PSymbols(fr, triv_vac = false)) : []
  pl = isempty(ps) ? [] : sort(collect(keys(ps)))
  pv = isempty(ps) ? [] : [ ps[l] for l in pl ]

  labels  = vcat(fl,rl,pl)
  nlabels = length(labels)

  zeropos    = [ x -> findfirst(x, labels) for x in zeros ]

  deleteat!(labels,zeropos)
  
  # functions to set up the rows of the gauge matrix
  sc = nzsc(fr)
  to_ind = Dict{Tuple{Int64,Int64,Int64},Int64}( sc[i] => i for i in 1:length(sc) ) 
  sr(x,y,z) = sparse_row(ZZ,[(to_ind((x,y,z)),1)])
  
  function gauge_row(l::FL)
    a,b,c,d,e = l[1:5]
    f = l[8]

    return sr(a,b,e) + sr(e,c,d) - sr(a,f,d) - sr(b,c,f)
  end
  
  function gauge_row(l::RL)
    a,b,c = l

    return sr(a,b,c) - sr(b,a,c)
  end

  d = conjugate_element(fr)
  function gauge_row(l::PL)
    a = l[1]
    
    return sr(d(a),a,1) - sr(a,d(a),1)
  end

  gaugemat = sparse_matrix(ZZ)

  for l in labels
    push!(gaugemat,gauge_row(l))
  end
  
  return gaugemat  
end

# adds columns containing to SMat until it has n columns
# adds no columns if n is smaller than or equal to
# number of columns of SMat.
# TODO: only works if SMat has ZZ elements
function pad!(m::SMat,n;kind="right")
  nc = number_of_columns(m) 
  nr = number_of_rows(m)
  diff = n - nc
  diff <= 0 && return m

  padmat = sparse_matrix(ZZ,nr,diff)

  if kind == "right"
    hcat(m,padmat)
  elseif kind == "left"
    # for some reason hcat does not work if you prepend a zero matrix...
    transpose(vcat(transpose(padmat),transpose(m))) 
  end
end

padleft(m::SMat,n)  = pad!(m,n,kind="left")
padright(m::SMat,n) = pad!(m,n,kind="right")

function d_sum(m1::SMat,m2::SMat)
  nc1 = number_of_columns(m1)
  nc2 = number_of_columns(m2)

  v1 = padright(m1,nc1+nc2)
  v2 = padleft(m2,nc1+nc2)
  vcat(v1,v2)
end

function gauge_split_basis( fr::FusionRing; zeros::Vector{FL} = [], includeonly::Vector{String} = ["FSymbols","PSymbols","RSymbols"])
  
end



#┏━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━┓
#┃                              to_unitary_gauge                                   ┃
#┗━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━┛

export to_unitary_gauge

"""to_unitary_gauge(fs::FSymbols)::Tuple{FSymbols,GaugeSymbols} returns
a tuple ( ufs, gs ) where ufs is an equivalent set of FSymbols for which the
F-matrices are unitary and gs is a set of gauge variables that perform
the transform from fs to ufs.
Returns ( nothing, nothing ) if a unitary gauge doesn't exist
"""
# Algorithm based on https://arxiv.org/pdf/2405.20075, pg 79
# We try to find a set of gauge symbols
function to_unitary_gauge(fs::FSymbols)

  return ifs = inverse_F_symbols(fs)
end

export to_unitary_gauge

function to_unitary_gauge(cat::SixJCategory)
    TCfs  = F_symbols(cat)
    field = base_ring(cat)
    emb   = cat.embedding


    # We will do everything over QQb and convert to a proper field
    # in the end

    # Note: we use guess so we have to check wheter the final result
    # is still a fusion category!

    deg = (degree ∘ minpoly ∘ gen ∘ base_ring)(cat)
    ϕ   = hom( field, QQb, guess( QQb, emb.r, deg ) )


    # convert TensorCategories F-Symbols to our F-symbols

    # convert TC label to 10-tuple (FLabel). Conventions
    # for 5th & 6th label differ between TC and Lyctr
    function label_to_tt(v::Vector{Int})::FL
        ( v[1], v[2], v[3], v[4], v[5], 1, 1, v[6], 1, 1 )
    end

    # embed the symbols and create a dict with FSL's as keys
    d = Dict( label_to_tt(k) => ϕ(v) for (k,v) in TCfs )

    mt = multiplication_table(cat)
    fs = F_symbols( d, QQb, mt; allvaluesknown = true )

    # invert the fsymbols
    ifs = inverse_F_symbols(fs)

    # set up a ring in squares of gauge symbols h(a,b,c) = g(a,b,c)^2
    fr = fusion_ring(mt)

    # nonzero-structure-constants
    sc = nzsc(fr)

    R, h = polynomial_ring( QQb, "𝒽" => sc )

    hdict = Dict( sc .=> h )
    𝒽(x,y,z) = hdict[(x,y,z)]

    pols = []


    invlabel(l::FL)::FL = l[[1, 2, 3, 4, 8, 9, 10, 5, 6, 7]]

    for l in labels(fs)
        a,b,c,d,f,_,_,e,_,_ = l
        push!(
            pols,
            𝒽(a,b,e) * 𝒽(e,c,d) * conj( value( fs, l ) ) -
            𝒽(a,f,d) * 𝒽(b,c,f) * value( ifs, invlabel(l) )
        )
    end

    sys = pol_sys( QQb, pols, constraint(R,[]), hdict )
end
