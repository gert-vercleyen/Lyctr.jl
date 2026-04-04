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

function pol_sys(
    ring, 
    pols, 
    constr, 
    knowns; 
    symm = missing, 
    attr = missing
    )
    PolSys(
        ring,
        pols,
        constr,
        knowns,
        symm,
        attr
    )
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

# TODO: implement ParametrizedPolSys for polynomial systems that can be iterated over but 
# are too large to store all at once. This would be useful for when we solve a subsystem using the smith decomposition
# and we get a large amount of solutions. 