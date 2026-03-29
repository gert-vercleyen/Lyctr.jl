

function is_inconsistent( s::PolSys )::Bool
  # Check whether pols contains nonzero pol
  # Check whether constr is not false
end

# TODO: The polynomials, constraints, and known values will be updated but for user defined attributes we should allow 
# users to add a way to update these as well 
"""evaluate( s::PolSys, d ) - Returns a polynomial system where the values of the variables have been updated via the dictionary d.
"""
function evaluate( s::PolSys, d )::PolSys

end


"""update_number_field( s::PolSys, ϕ ) - Returns a polynomial system where the definitions of the number fields have been updated 
to the be the image of ϕ and all values and coefficients belong to the image of ϕ.
The polynomials, constraints, and known values will be updated but user defined attributes will not.
"""
function update_number_field( s::PolSys, ϕ )::PolSys

end


