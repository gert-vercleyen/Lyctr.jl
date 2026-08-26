# Code for setting up equations arising from commutative diagrams
# We allow the equations to be set up over a ring for generality

# Pentagon equations 
export pentagon_system

# default ring is ℚ (ring should be expanded as we go)
function pentagon_system(fr::FusionRing; triv_vac = true, symbol = :ℱ)::PolSys
  return pentagon_system(fr, ℚ, QQemb; triv_vac = triv_vac, symbol = symbol)
end

function pentagon_system(fr::FusionRing, K::Field; triv_vac = true, symbol = :ℱ)::PolSys
  return pentagon_system(fr, K, nothing; triv_vac = triv_vac, symbol = symbol)
end

function pentagon_system(fr::FusionRing, K::Field, emb::CEmbOrNothing; triv_vac = true, symbol = :ℱ)::PolSys
  fs = F_symbols(fr; field = K, symbol = symbol, triv_vac = triv_vac)
  R  = polynomial_ring(fs)

  # construct pentagon polynomials
  pp = pent_polynomials(fs)

  # construct the invertibility constraints
  ic = invertibility_constraints(fs)

  # TODO: need to add snake cond [F^{aa*a}_a]^1_1 ≠ 0

  vrs = dict(fs)

  # set up the symmetries 
  s = gauge_symmetries(fs)

  fsymb(v::Vector{Int64}) = 1 ∈ v[1:3] ? R(1) : dict(v)

  return pol_sys(pp, ic, vrs; symm = s, emb = emb)
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

# TODO: don't distinguish between multiplicity and mf
# TODO: they give wrong result
# multiplicity-free pentagon equations
#

# convention F-symbols
#
# a   b   c                                        a   b   c
#  \ /   /                                          \   \ /
#   α   /                                            \   δ
#   e\ /      = \sum_{f,γ,δ} F[a,b,c,d,e,α,β,f,γ,δ]   \ /f
#     β                                                γ
#     |                                                |
#     d                                                d
#
function pent_polynomials(fsymb::FSymbols)::Vector{RingElem}
  fr     = fusion_ring(fsymb)
  r      = rank(fr)
  mt     = multiplication_table(fr)
  dic    = dict(fsymb)
  K      = polynomial_ring(fsymb)
  labels = keys(dic)

  # if trivial cat 
  r == 1 && return [dic[(1, 1, 1, 1, 1, 1, 1, 1, 1, 1)] - K(1)]

  fs(v::FL) = !haskey(dic, v) ? K(0) : dic[v]

  # Construct the polynomials
  pols = RingElem[]

  # add polynomials for which pent eqn has nonzero LHS
  for lab1 in labels
    f, c, d, e, g, _, _, l, _, _ = lab1

    function is_match(v::FL)
      return v[3] == l && v[4] == e && v[5] == f
    end

    matches = filter(is_match, labels)

    for lab2 in matches
      a, b, _, _, _, _, _, k, _, _ = lab2
      mα = mt[a,b,f]; mβ = mt[f,c,g];
      mγ = mt[g,d,e]; mε = mt[c,d,l];
      mθ = mt[a,k,e]; mη = mt[b,l,k];
      mδ = mt[f,l,e] # ζ in original paper but ζ is global gen of Qab
      for α ∈ 1:mα, β ∈ 1:mβ, γ ∈ 1:mγ, ε ∈ 1:mε, θ ∈ 1:mθ, η ∈ 1:mη
        mι(h) = mt[b,c,h]
        mκ(h) = mt[a,h,g]
        mλ(h) = mt[h,d,k]

        fac1(δδ)       = fs((f, c, d, e, g, β, γ, l, δδ, ε))
        fac2(δδ)       = fs((a, b, l, e, f, α, δδ, k, θ, η))
        fac3(hh,κκ,ιι) = fs((a, b, c, g, f, α, β, hh, κκ, ιι))
        fac4(hh,κκ,λλ) = fs((a, hh, d, e, g, κκ, γ, k, θ, λλ))
        fac5(hh,ιι,λλ) = fs((b, c, d, k, hh, ιι, λλ, l, η, ε))

        pol = sum( fac1(δ) * fac2(δ) for δ ∈ 1:mδ ) -
          sum(
            sum(
              fac3(h,κ,ι) *
              fac4(h,κ,λ) *
              fac5(h,ι,λ) for  ι ∈ 1:mι(h), κ ∈ 1:mκ(h), λ in 1:mλ(h);
                init = K(0)
                ) for h ∈ 1:r;
                  init = K(0)
          )

        if pol != K(0)
          push!(pols, pol)
        end
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
      mα = mt[a,b,f]; mβ = mt[f,c,g];
      mγ = mt[g,d,e]; mε = mt[c,d,l];
      mθ = mt[a,k,e]; mη = mt[b,l,k];
      mδ = mt[f,l,e]
      for α ∈ 1:mα, β ∈ 1:mβ, γ ∈ 1:mγ, ε ∈ 1:mε, θ ∈ 1:mθ, η ∈ 1:mη
        mι(h) = mt[b,c,h]
        mκ(h) = mt[a,h,g]
        mλ(h) = mt[h,d,k]

        fac3(hh,κκ,ιι) = fs((a, b, c, g, f, α, β, hh, κκ, ιι))
        fac4(hh,κκ,λλ) = fs((a, hh, d, e, g, κκ, γ, k, θ, λλ))
        fac5(hh,ιι,λλ) = fs((b, c, d, k, hh, ιι, λλ, l, η, ε))

        pol =
          sum(
            sum(
              fac3(h,κ,ι) *
              fac4(h,κ,λ) *
              fac5(h,ι,λ) for  ι ∈ 1:mι(h), κ ∈ 1:mκ(h), λ in 1:mλ(h);
                init = K(0)
                ) for h ∈ 1:r;
                  init = K(0)
          )

        pol != K(0) && push!(pols, pol)

      end
    end
  end

  return pols
end

function normalize_pol(pol)
  l = leading_coefficient(pol)
  if l != 1
    return pol/l
  else
    return pol
  end
end

function zero_struct_const(fr::FusionRing)::Vector{Tuple{Int64, Int64, Int64}}
  sc = nonzero_structure_constants(fr)
  r  = rank(fr)

  result = Tuple{Int64, Int64, Int64}[]

  for i in 1:r, j in 1:r, k in 1:r
    v = (i, j, k)
    if v ∉ sc
      push!(result, v)
    end
  end

  return result
end

export check_pentagon_equations

"""check_pentagon_equations(fs::FSymbols, rejectnonconstants=true)

Returns a couple ( b, pol_or_nothing ) where
* b is true and ind_or_nothing is nothing if the pentagon equations are satisfied
* b is false and pol_or_nothing is the first polynomial that is non-zero

By default a an equation is considered false if, after evaluation, its
polynomial is either non-constant (e.g. when not all values of the F-symbols
are known) or constant but non-zero.
By setting rejectnonconstants=false only non-zero constant polynomials will
lead to a value of b=false being returned.
"""

# sadly the quickest way is to check during construction so that we don't
# have to set up the whole system if the solution is incorrect
function check_pentagon_equations(fsymb::FSymbols, rejectnonconstants=true)
  fr     = fusion_ring(fsymb)
  r      = rank(fr)
  mt     = multiplication_table(fr)
  dic    = dict(fsymb)
  K      = polynomial_ring(fsymb)
  labels = collect(keys(dic))
  rnc    = rejectnonconstants
  formal_fs = Dict( v => k for (k,v) in invdict(fsymb) )


  fs(v::FL) = !haskey(dic, v) ? K(0) : formal_fs[v]

  # Construct the polynomials
  pols = RingElem[]

  # add polynomials for which pent eqn has nonzero LHS
  nl1 = length(labels)
  p_outer = Progress(nl1; desc = "1st F-symbol short path")
  for i in 1:nl1
    f, c, d, e, g, _, _, l, _, _ = labels[i]

    function is_match(v::FL)
      return v[3] == l && v[4] == e && v[5] == f
    end

    matches = filter(is_match, labels)
    nl2 = length(matches)
    
    p_inner = Progress(nl2; desc = "2nd F-symbol short path", offset = 1)
    for j in 1:nl2
      a, b, _, _, _, _, _, k, _, _ = matches[j]
      mα = mt[a,b,f]; mβ = mt[f,c,g];
      mγ = mt[g,d,e]; mε = mt[c,d,l];
      mθ = mt[a,k,e]; mη = mt[b,l,k];
      mδ = mt[f,l,e] # ζ in original paper but ζ is global gen of Qab
      for α ∈ 1:mα, β ∈ 1:mβ, γ ∈ 1:mγ, ε ∈ 1:mε, θ ∈ 1:mθ, η ∈ 1:mη
        mι(h) = mt[b,c,h]
        mκ(h) = mt[a,h,g]
        mλ(h) = mt[h,d,k]

        fac1(δδ)       = fs((f, c, d, e, g, β, γ, l, δδ, ε))
        fac2(δδ)       = fs((a, b, l, e, f, α, δδ, k, θ, η))
        fac3(hh,κκ,ιι) = fs((a, b, c, g, f, α, β, hh, κκ, ιι))
        fac4(hh,κκ,λλ) = fs((a, hh, d, e, g, κκ, γ, k, θ, λλ))
        fac5(hh,ιι,λλ) = fs((b, c, d, k, hh, ιι, λλ, l, η, ε))

        lhs = sum( fac1(δ) * fac2(δ) for δ ∈ 1:mδ ) 
        rhs = sum(
            sum(
              fac3(h,κ,ι) *
              fac4(h,κ,λ) *
              fac5(h,ι,λ) for  ι ∈ 1:mι(h), κ ∈ 1:mκ(h), λ in 1:mλ(h);
                init = K(0)
                ) for h ∈ 1:r;
                  init = K(0)
          )

        pol = lhs - rhs

        check = evaluate(pol,fsymb)
        if rnc
          if check != K(0)
            println("NONZERO LHS POL")
            @info "a, b, c, d, e, f, g, k, l: " ( a, b, c, d, e, f ,g, k, l )
            @info "lhs" evaluate( lhs, fsymb )
            @info "rhs" evaluate( rhs, fsymb )
            return (false,pol)
            end
        else
          is_constant(check) && check != K(0) && return (false,pol)
        end
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
      mα = mt[a,b,f]; mβ = mt[f,c,g];
      mγ = mt[g,d,e]; mε = mt[c,d,l];
      mθ = mt[a,k,e]; mη = mt[b,l,k];
      mδ = mt[f,l,e]
      for α ∈ 1:mα, β ∈ 1:mβ, γ ∈ 1:mγ, ε ∈ 1:mε, θ ∈ 1:mθ, η ∈ 1:mη
        mι(h) = mt[b,c,h]
        mκ(h) = mt[a,h,g]
        mλ(h) = mt[h,d,k]

        fac3(hh,κκ,ιι) = fs((a, b, c, g, f, α, β, hh, κκ, ιι))
        fac4(hh,κκ,λλ) = fs((a, hh, d, e, g, κκ, γ, k, θ, λλ))
        fac5(hh,ιι,λλ) = fs((b, c, d, k, hh, ιι, λλ, l, η, ε))

        pol =
          sum(
            sum(
              fac3(h,κ,ι) *
              fac4(h,κ,λ) *
              fac5(h,ι,λ) for  ι ∈ 1:mι(h), κ ∈ 1:mκ(h), λ in 1:mλ(h);
                init = K(0)
                ) for h ∈ 1:r;
                  init = K(0)
          )

        check = evaluate(pol,fsymb)
        if rnc
          if check != K(0)
            println("NONZERO LHS POL")
            @info "a, b, c, d, e, f, g, k, l: " ( a, b, c, d, e, f ,g, k, l )
            @info "rhs" evaluate( rhs, fsymb )
            return (false,pol)
            end
        else
          is_constant(check) && check != K(0) && return (false,pol)
        end
      end
    end
  end

  return ( true, nothing )
end

# TODO: add snake constraint or prove its unnecessary
function invertibility_constraints(fs::FSymbols)
  notzero(m) = inequality( det(m), ≠)
  fmats = F_matrices(fs)

  constr = filter(c -> !is_true(c), [notzero(m) for m in fmats])

  # snake constraint
  #snake(i) = inequality( (), ≠ )

  return constraint( polynomial_ring(fs), embedding(fs) , [[c] for c in constr])
end


export hexagon_system


function hexagon_system(fr::FusionRing; triv_vac = true, symbol = :ℛ, fsymbolvalues=nothing)::PolSys
  return pentagon_system(fr, ℚ, QQemb; triv_vac = triv_vac, symbol = symbol)
end

function hexagon_system(fr::FusionRing, K::Field; triv_vac = true, symbol = :ℛ, fsymbolvalues=nothing)::PolSys
  return pentagon_system(fr, K, nothing; triv_vac = triv_vac, symbol = symbol)
end


function hexagon_system(fr::FusionRing, K::Field, emb::CEmbOrNothing; triv_vac = true, symbols = (:ℱ,:ℛ), fsymbolvalues=nothing)::PolSys
  # TODO: need to take known and partially known F-values into account
  rs = R_symbols(fr; field = K, symbol = symbol, triv_vac = triv_vac)
  R  = polynomial_ring(rs)

  # construct pentagon polynomials
  hp = hex_polynomials(rs)

  # construct the invertibility constraints
  ic = invertibility_constraints(rs)

  vrs = dict(rs)

  # set up the symmetries 
  s = gauge_symmetries(rs)

  fsymb(v::Vector{Int64}) = 1 ∈ v[1:3] ? R(1) : dict(v)

  return pol_sys(pp, ic, vrs; symm = s)
end

# 
# multiplicity-free hexagon system
function hex_polynomials(rsymb::RSymbols,fsymb::FSymbols)::PolSys
  base_field(rsymb) != base_field(fsymb) && error("Number fields of symbols must match")

  symb = merge(fsymb,rsymb)
  
  R = polynomial_ring(symb)
  d = dict(symb)
  
  s(t::NTuple) = haskey(d,t) ? d(t) : R(0) 

  # Equation 1
  function fac11(L,G)
    rs((c,a,e,α,L)) *
    fs((a,c,b,d,e,L,β,g,μ,G)) *
    rs((c,b,g,G,ν))
  end

  function fac12(F,S,D,P)
    fs((c,a,b,d,e,α,β,F,S,D)) *
    rs((c,F,d,S,P)) *
    fs((a,b,c,d,F,D,P,g,μ,ν))
  end

  # Equation 2
  function fac21(L,G)
    irs((a,c,e,α,L)) *
    fs((a,c,b,d,e,L,β,g,μ,G)) *
    irs((b,c,g,G,ν))
  end
  
  function fac22(F,S,D,P)
    fs((c,a,b,d,e,α,β,F,S,D)) *
    irs((F,c,d,S,P)) *
    fs((a,b,c,d,F,D,P,g,μ,ν))
  end

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

function pivotal_equations(fr::FusionRing, k::Field)::PolSys
  return multiplicity(fr) === 1 ? mf_piv_eqns(fr, k) : piv_eqns(fr, k)
end

# multiplicity-free pivotal equations 
function mf_piv_eqns(fr::FusionRing, k::Field)::PolSys
  return print("Not implemented yet.")
  #TODO: implement
end

# pivotal equations with multiplicity
function piv_eqns(fr::FusionRing, k::Field)::PolSys
  return print("Not implemented yet.")
  #TODO: implement
end
