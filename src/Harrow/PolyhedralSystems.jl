# Code for setting up equations arising from commutative diagrams

# Pentagon equations 
export pentagon_equations

# default field is algebraic_closure of ℚ
function pentagon_equations(fr::FusionRing)::PolSys
    pentagon_equations( fr, ℚb )
end

function pentagon_equations(fr::FusionRing, k::Field )::PolSys
    multiplicity(fr) === 1 ? mf_pent_eqns(fr,k) : pent_eqns(fr,k) 
end

# multiplicity-free pentagon equations 
function mf_pent_eqns(fr::FusionRing,k::Field)::PolSys
    print("Not implemented yet.")

    #TODO: implement
end

# pentagon equations with multiplicity
function pent_eqns(fr::FusionRing,k::Field)::PolSys
    print("Not implemented yet.")
    #TODO: implement
end

export hexagon_equations

# default field is algebraic_closure of ℚ
function hexagon_equations(fr::FusionRing)::PolSys
    hexagon_equations( fr, ℚb )
end

function hexagon_equations(fr::FusionRing, k::Field )::PolSys
    multiplicity(fr) === 1 ? mf_hex_eqns(fr,k) : hex_eqns(fr,k) 
end

# multiplicity-free hexagon equations 
function mf_hex_eqns(fr::FusionRing,k::Field)::PolSys
    print("Not implemented yet.")

    #TODO: implement
end

# hexagon equations with multiplicity
function hex_eqns(fr::FusionRing,k::Field)::PolSys
    print("Not implemented yet.")

    #TODO: implement
end


export pivotal_equations

# default field is algebraic_closure of ℚ
function hexagon_equations(fr::FusionRing)::PolSys
    pivotal_equations( fr, ℚb )
end

function pivotal_equations(fr::FusionRing, k::Field )::PolSys
    multiplicity(fr) === 1 ? mf_piv_eqns(fr,k) : piv_eqns(fr,k) 
end

# multiplicity-free pivotal equations 
function mf_piv_eqns(fr::FusionRing,k::Field)::PolSys
    print("Not implemented yet.")
    #TODO: implement
end

# pivotal equations with multiplicity
function piv_eqns(fr::FusionRing,k::Field)::PolSys
    print("Not implemented yet.")
    #TODO: implement
end