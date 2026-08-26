#┏━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━┓
#┃                                  Abreviations                                   ┃
#┗━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━┛

CEmb = MapFromFunc{AbsSimpleNumField,QQBarField}

# If you want to work over abstract fields or fields unrelated to ℂ
CEmbOrNothing = Union{CEmb,NumFieldHom,Nothing}

#┏━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━┓
#┃                                   Inequality                                    ┃
#┗━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━┛


# Inequality of the form pol <s> 0 where s can be
# 1 : >,
# 2 : ≥,
# 3 : ≠
#
# wantembedded is necessary 
struct Inequality
  pol::RingElem      # Polynomial 
  emb::CEmbOrNothing # Embedding of base_ring(pol) into ℂ or nothing 
  t::Int64           # Type of inequality: 1 ↔ >, 2 ↔ ≥, 3 ↔ ≠
end

function inequality(polynomial::RingElem, type::Function; embedding=nothing)::Inequality
  return Inequality(polynomial, embedding, typetoint(type) )
end


function typetoint( type::Function )::Int64
  (type == >) && return 1
  (type == ≥) && return 2
  (type == ≠) && return 3
  error("Type must be either >, ≥, or ≠")
end

function polynomial(ineq::Inequality)::RingElem
  return ineq.pol
end

function type(ineq::Inequality)::Int64
  return ineq.t
end

function embedding(ineq::Inequality)::CEmbOrNothing
  return ineq.emb
end

function printstring(ineq::Inequality)
  t = type(ineq)
  s = if t == 1
    " > "
  elseif t == 2
    " ≥ "
  elseif t == 3
    " ≠ "
  end
  return string(polynomial(ineq)) * s * "0"
end

function Base.show(io::IO, ineq::Inequality)
  return print(io, printstring(ineq))
end

# TODO: should only be able to evaluate if
# 1. field of vals matches base_ring(polynomial(ineq))
# 2. field of vals extends base_ring(polynomial(ineq)) and
#    an embedding is given
function evaluate(
  ineq::Inequality, vars::Vector{Int}, vals::Vector{T}
)::Inequality where {T <: RingElement}
  return inequality( evaluate(polynomial(ineq), vars, vals), type(ineq))
end

function is_true(ineq::Inequality)
  p = polynomial(ineq)
  !is_constant(p) && return false

  t = type(ineq)

  t == 1 && return p > 0
  t == 2 && return p >= 0
  return p != 0
end

#┏━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━┓
#┃                                   Constraint                                    ┃
#┗━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━┛
# Constraint is a struct meant to provide a way to deal with a CNF form of
# propositions that take the form of inequalities.
# It solves the issue that there's no native OSCAR support to express and simplify
# propositions of the form
#   x[1] * x[3] - x[2] * x[4] ≠ 0 ∧ ( x[1] ≠ 0 ∨ x[2] ≠ 0)

# TODO: we probably don't need the polring since it should be stored in lis, even when lis is empty
struct Constraint
  polring::Ring                   # Polynomial ring of the variables appearing in the constraint
  emb::CEmbOrNothing              # Embedding of the base_ring(polring) into ℂ
  lis::Vector{Vector{Inequality}} # Vector of vector of inequalities
end

export constraint

function constraint(R::Ring, emb, list)::Constraint
  return Constraint(R, emb, list)
end

function constraint(R::Ring, list)::Constraint
  emb = ring == QQ ? QQemb : nothing
    
  return Constraint(R, emb, list)
end

export list

function list(c::Constraint)::Vector{Vector{Inequality}}
  return c.lis
end

export polynomial_ring

function polynomial_ring(c::Constraint)::Ring
  return c.polring
end

function embedding(c::Constraint)::CEmbOrNothing
  c.emb
end

export evaluate 

function evaluate(c::Constraint,vars,values)
  constraint
end
  

function Base.show(io::IO, constraint::Constraint)
  lis = list(constraint)

  if isempty(lis)
      println("true")
      return
  end

  function boolstring(b::Bool)
    return b ? "⊤" : "⊥"
  end

  function ineqstring(ineq)::String
    return typeof(ineq) === Bool ? boolstring(ineq) : printstring(ineq)
  end

  function orstring(l)::String
    return join(ineqstring.(l), " ∨ ")
  end


  orstrings = orstring.(lis)

  n = size(orstrings, 1)

  if n == 1
    print(orstrings[1])
    return
  end

  if n == 2
    print( orstrings[1] * " ∧ " * orstrings[2] )
    return
  end

  #for more elements we try to keep length under 80
  l1 = length(orstrings[begin])
  l2 = length(orstrings[end])

  startlength = n == 1 ? l1 : l1 + l2

  if startlength > 120
    str = join([orstrings[begin], "…", orstrings[end]], " ∧ ")
    print(io, str)
  else
    i = 2
    nextstring = orstrings[i]
    l = startlength + length(nextstring)

    strings = [orstrings[begin]]
    while l < 120 && i < n
      push!(strings, orstrings[i])
      i = i + 1
      l = l + length(orstrings[i])
    end

    if size(strings, 1) < n - 1
      push!(strings, "…")
      push!(strings, orstrings[end])
    end

    if size(strings,1) == n - 1
      push!(strings, orstrings[end])
    end

    print(io, join(strings, " ∧ "))
  end
end

#┏━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━┓
#┃                                     PolSys                                      ┃
#┗━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━┛

# This struct should be as optimal as possible. It might be best to
# parametrize the number field so we can give predefined types to all struct fields
struct PolSys
  pls::Vector{<:RingElem} # Polynomials of the system. Could be Laurent polynomials
  emb::CEmbOrNothing      # Embedding of the base_ring(polring) into ℂ or nothing
  cnstr::Constraint       # Constraints of the system
  knwns::Dict{Any,Any}    # Dictionary of known values. 
  symm                    # Symmetries of the variables
  __attr                  # Further info
end

export pol_sys

function pol_sys(pols, constr, knowns; symm = missing, attr = missing, emb = nothing)
  isempty(knowns) && error("The dictionary of variables should be non-empty")
  return PolSys(pols, emb, constr, knowns, symm, attr)
end

export polynomials

function polynomials(s::PolSys)
  return s.pls
end

export constraint

function constraint(s::PolSys)
  return s.cnstr
end

export known_values

function known_values(s::PolSys)
  return s.knwns
end

export polynomial_ring

function polynomial_ring(s::PolSys)
  return (parent ∘ first ∘ values ∘ known_values)(s)
end

export base_field

function base_field(s::PolSys)
  return base_ring(polynomial_ring(s))
end


function embedding(s::PolSys)
  return s.emb
end

# TODO: implement ParametrizedPolSys for polynomial systems that can be iterated 
# over but are too large to store all at once. This would be useful for when we 
# solve a subsystem using the smith decomposition and we get a large amount of 
# solutions. 

function Base.show(io::IO, ps::PolSys)
  pols = polynomials(ps)
  lps  = string(length(pols))
  cns  = constraint(ps)
  lcns = string(length(list(cns)))
  vrs  = known_values(ps)
  lvrs = string(length(vrs))
  fld  = string(base_field(ps))
  emb  = string(embedding(ps))

  function showpolies(polies)
    if length(polies) < 8
      println.(polies)
    else
      println.(polies[1:3])
      println("⋮")
      println.(polies[end-2:end])
    end
  end

  function showvars(vars)
    l = sort(collect(vars), by = first)
    function printpair(p)
      k, v = p
      print(k)
      print(" => ")
      println(v)
    end

    if length(l) < 8
      printpair.(l)
    else
      printpair.(l[1:3])
      println("⋮")
      printpair.(l[end-2:end])
    end
  end

  str = "Collection of " * lps * " polynomials and " * lcns * " inequalities in " * lvrs * " variables over " * fld

  println(str)
  if emb != "nothing"
    println("embedded via \n", emb)
  else
    println("without embedding.")
  end
  println("===Polynomials===")
  showpolies(pols)
  println("\n===Constraints===")
  show(cns)
  println("\n\n===Variables===")
  showvars(vrs)
end
