# TODO:
# This struct should be as optimal as possible. It might be best to
# parametrize the number field so we can give predefined types to all struct fields
struct PolSys
  rng::Ring            # base ring over which pols are defined
  pls            # Polynomials of the system
  cnstr          # Constraints of the system
  knwns          # Dictionary of known values
  symm           # Symmetries of the variables
  __attr         # Further info
end

function pol_sys(ring, pols, constr, knowns; symm = missing, attr = missing)
  return PolSys(ring, pols, constr, knowns, symm, attr)
end

function base(s::PolSys)
  return s.rng
end

function polynomials(s::PolSys)
  return s.pls
end

function constraints(s::PolSys)
  return s.cnstr
end

function known_values(s::PolSys)
  return s.knwns
end

# TODO: implement ParametrizedPolSys for polynomial systems that can be iterated over but 
# are too large to store all at once. This would be useful for when we solve a subsystem using the smith decomposition
# and we get a large amount of solutions. 

# Inequality of the form pol <s> 0 where s can be
# 1 : >, 
# 2 : ≥, 
# 3 : ≠

struct Inequality
  rng::Ring
  pol
  t::Int64
end

function inequality( ring, polynomial, type )
  s = if type == >
    1
  elseif type == ≥
    2
  elseif type == ≠
    3
  end

  Inequality( ring, polynomial, s )
end

polynomial( ineq::Inequality ) = ineq.pol

type( ineq::Inequality ) = ineq.t

ring( ineq::Inequality ) = ineq.rng

function Base.show( io::IO, ineq::Inequality )
  t = type(ineq)
  s = if t == 1
    " > "
  elseif t == 2
    " ≥ "
  elseif t == 3
    " ≠ "
  end
  print( io, string(polynomial(ineq)) * s * "0" )
end

function evaluate( ineq::Inequality, vars::Vector{Int}, vals::Vector{T} ) where {T <: RingElement}
  inequality(
    ring(ineq),
    evaluate( polynomial(ineq), vars, vals ),
    type(ineq)
  )
end

function is_true( ineq::Inequality ) 
  p = polynomial(ineq)
  !is_constant(p) && return false
  
  t = type(ineq)
  if t == 1
    return p > 0 
  elseif t == 2
    return p >= 0
  else 
    return p != 0
  end 
end