# Code for setting up equations arising from commutative diagrams
# We allow the equations to be set up over a ring for generality

# Pentagon equations 
export pentagon_system

# default ring is ℚ (ring should be expanded as we go)
function pentagon_system(fr::FusionRing; triv_vac = true, symbol = :ℱ)::PolSys
  return pentagon_system(fr, ℚ; triv_vac = triv_vac, symbol = symbol)
end

function pentagon_system(fr::FusionRing, k::Ring; triv_vac = true, symbol = :ℱ)::PolSys
  if multiplicity(fr) === 1
    mf_pent_sys(fr, k; triv_vac = triv_vac, symbol = :ℱ)
  else
    pent_sys(fr, k; triv_vac = triv_vac, symbol = :ℱ)
  end
end

function mf_pent_sys(fr::FusionRing, K::Ring; symbol = :ℱ, triv_vac = true)::PolSys
  r = rank(fr)
  R, dict = F_symbols(fr; ring = K, symbol = symbol)

  # construct pentagon polynomials
  pp = mf_pent_polynomials(fr, dict, K; triv_vac = triv_vac)

  # construct the invertibility constraints
  ic = mf_invertibility_constraints(R, dict; triv_vac = triv_vac)

  # set up the symmetries 
  s = mf_gauge_symmetries(r, dict, K; triv_vac = triv_vac)

  fsymb(v::Vector{Int64}) = 1 ∈ v[1:3] ? R(1) : dict(v)

  if triv_vac
    set_triv_vac(k,v) = 1 ∈ k[1:3] ? 1 : v 

    variables = Dict(k => set_triv_vac(k,v) for (k,v) in dict )

    return pol_sys(R, pp, ic, variables; symm = s)
  else
    variables = dict 

    return pol_sys(R, pp, ic, variables; symm = s)
  end

end

# TODO: check whether these are actually sparse. 
# I have the feeling they might not be due to how 
# slow the algorithm is 

#= TODO: check whether following construction is faster
julia> R, (x, y) = polynomial_ring(ZZ, [:x, :y])
(Multivariate polynomial ring in 2 variables over integers, AbstractAlgebra.Generic.MPoly{BigInt}[x, y])

julia> C = MPolyBuildCtx(R)
Builder for an element of R

julia> push_term!(C, ZZ(3), [1, 2]);


julia> push_term!(C, ZZ(2), [1, 1]);


julia> push_term!(C, ZZ(4), [0, 0]);


julia> f = finish(C)
3*x*y^2 + 2*x*y + 4

julia> push_term!(C, ZZ(4), [1, 1]);


julia> f = finish(C)
4*x*y

julia> S, (x, y) = polynomial_ring(QQ, [:x, :y])
(Multivariate polynomial ring in 2 variables over rationals, AbstractAlgebra.Generic.MPoly{Rational{BigInt}}[x, y])

julia> f = S(Rational{BigInt}[2, 3, 1], [[3, 2], [1, 0], [0, 1]])
2*x^3*y^2 + 3*x + y
=#

# multiplicity-free pentagon equations 
function mf_pent_polynomials(fr, dict, K::Ring; triv_vac = true)
  r = rank(fr)

  # if trivial cat 
  if r == 1
    return [ dict[ [1,1,1,1,1,1,1,1,1,1] ] - K(1) ] 
  end

  labels = keys(dict)

  function fs(v::Vector{Int64})
    if triv_vac
      !haskey(dict, v) && return K(0)
      1 ∈ v[1:3] && return K(1)
      return dict[v]
    else
      return !haskey(dict, v) ? K(0) : dict[v]
    end
  end

  # Construct the polynomials

  pols = RingElem[]

  # add polynomials for which pent eqn has nonzero LHS
  for lab1 in labels
    f, c, d, e, g,_ ,_ , l, _, _ = lab1

    function is_match(v::Vector{Int64})
      return v[3] == l && v[4] == e && v[5] == f
    end

    matches = filter(is_match, labels)

    for lab2 in matches
      a, b, _, _, _, _, _, k, _, _ = lab2
      pol =
        fs([f, c, d, e, g, 1, 1, l, 1, 1 ]) * 
        fs([a, b, l, e, f, 1, 1, k, 1, 1 ]) - 
        sum(
          fs([a, b, c, g, f, 1, 1, h, 1, 1 ])*
          fs([a, h, d, e, g, 1, 1, k, 1, 1 ])*
          fs([b, c, d, k, h, 1, 1, l, 1, 1 ]) for
          h in 1:r
        )

      if pol != K(0)
        push!(pols, pol)
      end
    end
  end

  # add polynomials for which pent eqn has zero LHS and nonzero RHS
  # This is done by constructing the symmetric tree with non-existent 
  # bottom fusion channel N[f,l,e] and matching the other labels

  sc = nonzero_structure_constants(fr)
  zsc = zero_struct_const(fr)

  for n1 in zsc
    f, l, e = n1

    function is_match2(v)::Bool
      return v[3] == f
    end

    function is_match3(v)::Bool
      return v[3] == l
    end

    matches2 = filter(is_match2, sc)
    matches3 = filter(is_match3, sc)

    for n2 in matches2, n3 in matches3, k in 1:r, g in 1:r
      a, b, = n2
      c, d, = n3
      pol = sum(
        fs([a, b, c, g, f, 1, 1, h, 1, 1 ])*
        fs([a, h, d, e, g, 1, 1, k, 1, 1 ])*
        fs([b, c, d, k, h, 1, 1, l, 1, 1 ]) for
        h in 1:r
      )

      if pol != K(0)
        push!(pols, pol)
      end
    end
  end

  return pols
end

function normalize_pol( pol ) 
  l = leading_coefficient(pol)
  if l != 1
    return pol/l
  else 
    return pol
  end
end

function zero_struct_const(fr::FusionRing)::Vector{Tuple{Int64,Int64,Int64}}
  sc = nonzero_structure_constants(fr)
  r  = rank(fr)

  result = Tuple{Int64,Int64,Int64}[]

  for i in 1:r, j in 1:r, k in 1:r
    v = (i, j, k)
    if v ∉ sc
      push!(result, v)
    end
  end

  return result
end



function mf_invertibility_constraints(R::Ring, dict; triv_vac = triv_vac) 
  K = base_ring(R)
  notzero(m) = inequality( K, det(m), ≠ )
  _, _, fmats = F_matrices( R, dict, triv_vac = triv_vac )

  constr = filter( c -> !is_true(c),  [ notzero(m) for m in fmats ] )

  constraint( R, [ [ c ] for c in constr ]  )
end

function mf_gauge_symmetries(r, dict, K; triv_vac = triv_vac) end

# pentagon equations with multiplicity
function pent_eqns(fr::FusionRing, k::Ring)::PolSys
  return print("Not implemented yet.")
  #TODO: implement
end

export hexagon_equations

# default ring is algebraic_closure of ℚ
function hexagon_equations(fr::FusionRing)::PolSys
  return hexagon_equations(fr, ℚb)
end

function hexagon_equations(fr::FusionRing, k::Ring)::PolSys
  return multiplicity(fr) === 1 ? mf_hex_eqns(fr, k) : hex_eqns(fr, k)
end

# multiplicity-free hexagon equations 
function mf_hex_eqns(fr::FusionRing, k::Ring)::PolSys
  return print("Not implemented yet.")

  #TODO: implement
end

# hexagon equations with multiplicity
function hex_eqns(fr::FusionRing, k::Ring)::PolSys
  return print("Not implemented yet.")

  #TODO: implement
end

export pivotal_equations

# default ring is algebraic_closure of ℚ
function pivotal_equations(fr::FusionRing)::PolSys
  return pivotal_equations(fr, ℚb)
end

function pivotal_equations(fr::FusionRing, k::Ring)::PolSys
  return multiplicity(fr) === 1 ? mf_piv_eqns(fr, k) : piv_eqns(fr, k)
end

# multiplicity-free pivotal equations 
function mf_piv_eqns(fr::FusionRing, k::Ring)::PolSys
  return print("Not implemented yet.")
  #TODO: implement
end

# pivotal equations with multiplicity
function piv_eqns(fr::FusionRing, k::Ring)::PolSys
  return print("Not implemented yet.")
  #TODO: implement
end
