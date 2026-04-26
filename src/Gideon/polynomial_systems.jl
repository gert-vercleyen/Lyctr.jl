# Methods for manipulating polynomial systems 

# Checks whether a polynomial system is inconsistent in a trivial way
function is_trivially_inconsistent(s::PolSys)::Bool
  # Check whether pols contains nonzero pol
  is_nz_const_pol( p ) = is_constant(p) && constant_coefficient(p) != 0
  any( is_nz_const_pol, polynomials(s) ) && return true

  # Check whether constr is not false
  is_false_prop( p ) = p[1] === false
  any( is_false_prop, list(constraint(s)) ) && return true 

  return false
end

iti = is_trivially_inconsistent

function looks_fine(s::PolSys)::Bool
  return !is_trivially_inconsistent(s)
end

# TODO: The polynomials, constraints, and known values will be updated but for user defined attributes we should allow 
# users to add a way to update these as well 
"""evaluate( s::PolSys, d ) - Returns a polynomial system where the values of the variables have been updated via the dictionary d.
"""
function evaluate(s::PolSys, dict)::PolSys
  is_updated( ( k, v ) ) = k != v

  actualupdates = filter( is_updated, dict )

  # If no update required
  is_empty(actualupdates) && return s

  vars = collect(keys(actualupdates))
  vals = collect(values(actualupdates))

  evaluate( s, vars, vals )
end

function evaluate(s::PolSys, vars, vals )
  function update(p)
    return evaluate(p, vals, vars )
  end

  function mapupdate(list)
    return map( update, list )
  end

  newpols   = map( update, polynomials(s) )
  newconstr = constraint(R, map( mapupdate, list(constraint(s)) ) )
  newknowns = Dict( k => update(v) for (k,v) in known_values(s) )

  # TODO: implement updating of symmetries
  pol_sys( 
    polynomial_ring(s),
    newpols,
    newconstr,
    newknowns  
  )
end

"""update_number_field( s::PolSys, ϕ ) - Returns a polynomial system where the definitions of the number fields have been updated 
to the be the image of ϕ and all values and coefficients belong to the image of ϕ.
"""
function update_number_field(s::PolSys, ϕ)::PolSys
end