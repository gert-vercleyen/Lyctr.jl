# TODO:
# This struct should be as optimal as possible. It might be best to
# parametrize the number field so we can give predefined types to all struct fields

struct PolSys
       rng            # base ring over which pols are defined
       pls            # Polynomials of the system
       cnstr          # Constraints of the system
       knwns          # Dictionary of known values
       symm           # Symmetries of the variables
       __attr         # Further info
end

function base( s::PolSys )
    return s.rng
end

function polynomials( s::PolSys )
    return s.pls
end

function constraints( s::PolSys )
    return s.cnstr
end

function known_values( s::PolSys )
    return s.knwns
end

function is_inconsistent( s::PolSys )::Bool
  # Check whether pols contains nonzero pol
  # Check whether constr is not false
end

"""evaluate( s::PolSys, d ) - Returns a polynomial system where the values of the variables have been updated via the dictionary d.
The polynomials, constraints, and known values will be updated but user defined attributes will not.
"""
function evaluate( s::PolSys, d )::PolSys

end

function update_number_field( s::PolSys )::PolSys
end
