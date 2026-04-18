# Code for setting up equations arising from commutative diagrams
# We allow the equations to be set up over a ring for generality

# Pentagon equations 
export pentagon_system

# default ring is algebraic_closure of ℚ
function pentagon_system(fr::FusionRing; triv_vac = true, symbol = :ℱ)::PolSys
  return pentagon_system(fr, ℚb; triv_vac = triv_vac, symbol = symbol)
end

function pentagon_system(fr::FusionRing, k::Ring; triv_vac = true, symbol = :ℱ)::PolSys
  if multiplicity(fr) === 1
    mf_pent_sys(fr, k; triv_vac = triv_vac, symbol = :ℱ)
  else
    pent_sys(fr, k; triv_vac = triv_vac, symbol = :ℱ)
  end
end

function mf_pent_sys(fr::FusionRing, K::Ring; symbol = :ℱ, triv_vac = true)::PolSys
  # TODO: should only construct F-symbols and ring once
  r = rank(fr)
  R, dict = F_symbols(fr; ring = K, symbol = symbol)

  # construct pentagon polynomials
  pp = mf_pent_eqns(r, dict, K; triv_vac = triv_vac)

  # construct the invertibility constraints
  ic = mf_invertibility_constraints(r, dict, K; triv_vac = triv_vac)

  # set up the symmetries 
  s = mf_gauge_symmetries(r, dict, K; triv_vac = triv_vac)

  fsymb(v::Vector{Int64}) = 1 ∈ v[1:3] ? K(1) : dict(v)

  variables = Dict(key => fsymb(key) for key in keys(dict))

  return pol_sys(R, pp, ic, variables; symm = s)
end

# multiplicity-free pentagon equations 
function mf_pent_eqns(r::Int64, dict, K::Ring; triv_vac = true)

  # trivial cat 
  if r == 1
    return [ dict[ [1,1,1,1,1,1,1,1,1,1] ] - K(1) ] 
  end


  labels = keys(dict)

  function fs(v::Vector{Int64})
  if triv_vac
      1 ∈ v[1:3] && return K(1)
      return !haskey(dict, v) ? K(0) : dict[v]
    else
      return !haskey(dict, v) ? K(0) : dict[v]
    end
  end

  # Construct the polynomials

  pols = RingElem[]

  # add polynomials for which pent eqn has nonzero LHS
  for lab1 in labels
    f, c, d, e, g, l = lab1

    function is_match(v::Vector{Int64})
      return v[3] == l && v[4] == e && v[5] == f
    end

    matches = filter(is_match, labels)

    for lab2 in matches
      a, b, _, _, _, k = lab2
      pol =
        fs([f, c, d, e, g, l]) * fs([a, b, l, e, f, k]) - sum(
          fs([a, b, c, g, f, h])*fs([a, h, d, e, g, k])*fs([b, c, d, k, h, l]) for
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

  sc = non_zero_structure_constants(fr)
  zsc = zero_struct_const(fr)

  for n1 in zsc
    f, l, e = n1

    function is_match2(v::Vector{Int64})::Bool
      return v[3] == f
    end

    function is_match3(v::Vector{Int64})::Bool
      return v[3] == l
    end

    matches2 = filter(is_match1, sc)
    matches3 = filter(is_match2, sc)

    for n2 in matches2, n3 in matches3, k in 1:r, g in 1:r
      a, b, = n2
      c, d, = n3
      pol = sum(
        fsymb([a, b, c, g, f, h])*fsymb([a, h, d, e, g, k])*fsymb([b, c, d, k, h, l]) for
        h in 1:r
      )

      if pol != K(0)
        push!(pols, pol)
      end
    end
  end

  return pols
end

function zero_struct_const(fr::FusionRing)
  sc = non_zero_structure_constants(fr)
  r  = rank(fr)

  result = Vector{Int64}[]

  for i in 1:r, j in 1:r, k in 1:r
    v = [i, j, k]
    if v ∉ sc
      push!(result, v)
    end
  end

  return result
end

function mf_invertibility_constraints(r, dict, K; triv_vac = triv_vac) end

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
function hexagon_equations(fr::FusionRing)::PolSys
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
