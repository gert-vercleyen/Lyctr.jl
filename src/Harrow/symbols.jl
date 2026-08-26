# This file contains constructors for all kinds of symbols whose labels represent objects in a fusion category

#┏━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━┓
#┃                                    F_Labels                                     ┃
#┗━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━┛

export F_labels
# F_labels returns a list of 10-element lists corresponding to well-formed F-symbols
function F_labels(fr::FusionRing)::Vector{FL}
  r  = rank(fr)
  mt = multiplication_table(fr)
  sc = nzsc(fr)

  lbls = NTuple{10, Int64}[]

  for lab1 in sc, lab2 in sc
    if lab1[3] == lab2[1]
      a, b, e = lab1;
      m1      = mt[a, b, e]
      c, d    = lab2[2:end];
      m2      = mt[e, c, d]
      for f in 1:r
        m3 = mt[b, c, f]
        m4 = mt[a, f, d]
        if m3 * m4 != 0
          for α in 1:m1, β in 1:m2, γ in 1:m3, δ in 1:m4
            push!(lbls, (a, b, c, d, e, α, β, f, γ, δ))
          end
        end
      end
    end
  end

  return sort(lbls)
end

#┏━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━┓
#┃                                   F_symbols                                     ┃
#┗━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━┛
#
export F_symbols

"""F_symbols( fr::FusionRing; field=ℚ, embedding = QQemb symbol = :ℱ ) returns an object of type FSymbols  
"""
function F_symbols(fr::FusionRing; field = ℚ, embedding = QQemb, symbol = :ℱ, triv_vac = true)::FSymbols
  fl = F_labels(fr)
  n  = length(fl)

  R, f = polynomial_ring(field, symbol => 1:n; cached = true)

  function value(i::Int64)
    if triv_vac
      1 ∈ fl[i][1:3] ? R(1) : f[i]
    else
      f[i]
    end
  end

  dict  = Dict(fl[i] => value(i) for i in 1:n)
  idict = Dict( [ f[i] => fl[i] for i ∈ 1:n ]... )

  return F_symbols(dict, R, idict, embedding, multiplication_table(fr))
end

export vac_F_symbols

function vac_F_symbols(R::Field, dict)
  return (R, [v for (k, v) in dict if 1 ∈ k[1:3]])
end

export F_matrices

function F_matrices(
  fr::FusionRing; field = ℚ, embedding = QQemb, symbol = :ℱ, triv_vac = true
)
  fs = F_symbols(fr; field = field, embedding=embedding, symbol = symbol)
  return F_matrices(fs)
end

function F_matrices(fs::FSymbols)
  fr = fusion_ring(fs)

  ffmats = formal_F_matrices(F_labels(fr))

  return map( mat -> Feval(fs,mat), ffmats)
end

export formal_F_matrices

function formal_F_matrices(fr::FusionRing)::Vector{Matrix{FL}}
  return formal_F_matrices(F_labels(fr))
end

function formal_F_matrices(flabels::Vector{FL})::Vector{Matrix{FL}}
  grouped_labels = gather_by(l -> l[1:4], flabels)

  function formal_fmat(lbls::Vector{FL})::Matrix{FL}
    # associativity of fusion ring makes size of lbls square
    m = Int64(sqrt(size(lbls, 1)))
    return reshape(lbls, m, m)
  end

  return formal_fmat.(grouped_labels)
end


function Feval(fs::FSymbols, mat::Matrix{FL})
  return matrix(polynomial_ring(fs), value(fs,mat))
end

"""inverse_F_symbols(fs::FSymbols) returns the inverse F-symbols corresponding
to the F-symbols fs. The values of the inverse F-symbols belong to the fraction
field of the polynomial ring of the F-symbols. This function also works for F-symbols
whose values are not yet determined, in which case the formula for the inverse
of a general matrix is used. 
"""
function inverse_F_symbols(fs::FSymbols)
  # matrices of labels of F-symbols
  lmats = formal_F_matrices(fusion_ring(fs))
  # matrices of values of F-symbols
  R  = polynomial_ring(fs)
  FF = fraction_field(R)
  Fmats = map( m -> matrix( FF, value(fs, m) ), lmats )

  # since the labels of an F-matrix denote the anyon
  # types the label sets need to be swapped for inverse
  # F symbols to make sense
  swap_labels(tt::FL) = tt[[1, 2, 3, 4, 8, 9, 10, 5, 6, 7]]
  function swap_labels(m::Matrix{FL})
      [ swap_labels( m[i,j] ) for i in 1:nrows(m), j in 1:ncols(m) ]
  end

  invlmats = swap_labels.(lmats)

  invFmats = inv.(Fmats)

  d = Dict{FL,FieldElem}()
  for ind ∈ eachindex(invlmats)
    for i ∈ 1:nrows(invlmats[ind]), j ∈ 1:ncols(invlmats[ind])
      push!(d, invlmats[ind][i,j] => invFmats[ind][i,j] )
    end
  end
  
  return inverse_F_symbols(
    d,
    R,
    invdict(fs),
    multiplication_table(fs),
    allvaluesknown = allvaluesknown(fs),
    inunitarygauge = inunitarygauge(fs),
    inminfield     = inminfield(fs)
  )
end

#┏━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━┓
#┃                                    R_Labels                                     ┃
#┗━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━┛

export R_labels
# F_labels returns a list of 10-element lists corresponding to well-formed F-symbols
function R_labels(fr::FusionRing)::Vector{RL}
  r  = rank(fr)
  mt = multiplication_table(fr)
  sc = nzsc(fr)

  lbls = RL[]

  for lab1 in sc
    a, b, c = lab1;
    m       = mt[a, b, c]
    for α in 1:m, β in 1:m
      push!(lbls, (a, b, c, α, β ) )
    end
  end

  return sort(lbls)
end
#┏━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━┓
#┃                                   R_symbols                                     ┃
#┗━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━┛
#
export R_symbols

"""R_symbols( fr::FusionRing; field=ℚ, embedding = QQemb symbol = :ℛ ) returns an object of type RSymbols  
"""
function R_symbols(fr::FusionRing; field = ℚ, embedding = QQemb, symbol = :ℛ, triv_vac = true)::RSymbols
  rl = R_labels(fr)
  n  = length(rl)

  R, r = polynomial_ring(field, symbol => 1:n; cached = true)

  function value(i::Int64)
    if triv_vac
      1 ∈ rl[i][1:2] ? R(1) : r[i]
    else
      r[i]
    end
  end

  dict  = Dict( rl[i] => value(i) for i in 1:n)
  idict = Dict( [ r[i] => rl[i] for i ∈ 1:n ]... )

  return R_symbols(dict, R, idict, embedding, multiplication_table(fr))
end

export R_matrices

function R_matrices(
  fr::FusionRing; field = ℚ, embedding = QQemb, symbol = :ℛ, triv_vac = true
)
  rs = R_symbols(fr; field = field, embedding=embedding, symbol = symbol)
  return R_matrices(rs)
end

function R_matrices(rs::RSymbols)
  fr = fusion_ring(rs)

  frmats = formal_R_matrices(R_labels(fr))

  return map( mat -> Reval(rs,mat), frmats)
end

export formal_R_matrices

function formal_R_matrices(fr::FusionRing)::Vector{Matrix{RL}}
  return formal_R_matrices(R_labels(fr))
end

function formal_R_matrices(rlabels::Vector{RL})::Vector{Matrix{RL}}
  grouped_labels = gather_by(l -> l[1:3], rlabels)

  function formal_rmat(lbls::Vector{RL})::Matrix{RL}
    m = Int64(sqrt(size(lbls, 1)))
    return reshape(lbls, m, m)
  end

  return formal_rmat.(grouped_labels)
end


function Reval(rs::RSymbols, mat::Matrix{RL})
  return matrix(polynomial_ring(rs), value(rs,mat))
end

"""inverse_F_symbols(fs::FSymbols) returns the inverse F-symbols corresponding
to the F-symbols fs. The values of the inverse F-symbols belong to the fraction
field of the polynomial ring of the F-symbols. This function also works for F-symbols
whose values are not yet determined, in which case the formula for the inverse
of a general matrix is used. 
"""
function inverse_R_symbols(rs::FSymbols)
  # matrices of labels of R-symbols
  lmats = formal_R_matrices(fusion_ring(rs))
  # matrices of values of R-symbols
  R  = polynomial_ring(fs)
  FF = fraction_field(R)
  Rmats = map( m -> matrix( FF, value(rs, m) ), lmats )

  # since the labels of an R-matrix denote the anyon
  # types, the top two labels need to be swapped for inverse
  # R-symbols to make sense
  swap_labels(tt::FL) = tt[[2, 1, 3, 4, 5]]
  function swap_labels(m::Matrix{RL})
      [ swap_labels( m[i,j] ) for i in 1:nrows(m), j in 1:ncols(m) ]
  end

  invlmats = swap_labels.(lmats)

  invRmats = inv.(Rmats)

  d = Dict{RL,FieldElem}()
  for ind ∈ eachindex(invlmats)
    for i ∈ 1:nrows(invlmats[ind]), j ∈ 1:ncols(invlmats[ind])
      push!(d, invlmats[ind][i,j] => invRmats[ind][i,j] )
    end
  end
  
  return inverse_R_symbols(
    d,
    R,
    invdict(rs),
    multiplication_table(rs),
    allvaluesknown = allvaluesknown(rs),
    inunitarygauge = inunitarygauge(rs),
    inminfield     = inminfield(rs)
  )
end


#┏━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━┓
#┃                                   P_symbols                                     ┃
#┗━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━┛

export P_symbols

"""P_symbols( fr::FusionRing; field=ℚ, embedding = QQemb symbol = :𝒫 ) returns an object of type PSymbols  
"""
function P_symbols(fr::FusionRing; field = ℚ, embedding = QQemb, symbol = :𝒫, triv_vac = true)::PSymbols
  rk = rank(fr)
  pl = [ (i,) for i in 1:rk ]

  R, p = polynomial_ring(field, symbol => 1:rk; cached = true)

  function value(i::Int64)
    if triv_vac
      i == 1 ? R(1) : p[i]
    else
      p[i]
    end
  end

  dict  = Dict( pl[i] => value(i) for i in 1:rk)
  idict = Dict( [ p[i] => pl[i] for i ∈ 1:rk ]... )

  return P_symbols(dict, R, idict, embedding, multiplication_table(fr))
end

export P_matrices

function P_matrices(
  fr::FusionRing; field = ℚ, embedding = QQemb, symbol = :𝒫, triv_vac = true
)
  ps = P_symbols(fr; field = field, embedding=embedding, symbol = symbol)
  return P_matrices(ps)
end

function P_matrices(ps::PSymbols)
  frmats = formal_p_matrices(fusion_ring(ps))

  return map( mat -> Peval(ps,mat), frmats)
end

export formal_P_matrices

function formal_P_matrices(fr::FusionRing)::Vector{Matrix{PL}}
  return [ [ (i) ] for i in 1:rank(fr) ] 
end

function Peval(ps::PSymbols, mat::Matrix{PL})
  return matrix(polynomial_ring(ps), value(ps,mat))
end


#┏━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━┓
#┃                                    G-symbols                                    ┃
#┗━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━┛

# G-symbols are Lyctr's symbols used to express gauge degrees of freedom

export G_symbols

function G_symbols(fr::FusionRing; field = ℚ, embedding = QQemb, symbol = :𝒢)::GSymbols
  # labels of G symbols
  gl = GL[]
  r  = rank(fr)
  mt = multiplication_table(fr)
  for i ∈ 1:r, j ∈ 1:r, k ∈ 1:r
    mult = mt[i,j,k]
    for α ∈ 1:mult, β ∈ 1:mult
      push!(gl, (i,j,k,α,β))
    end
  end
  
  n = length(gl)

  R, g = polynomial_ring( field, symbol => 1:n, cached = true )
  
  dict = Dict( gl[i] => g[i] for i ∈ 1:n )

  idict = Dict( [ v => k for (k,v) ∈ dict ] ... )

  return G_symbols(dict,R,idict,embedding,mt)
end
