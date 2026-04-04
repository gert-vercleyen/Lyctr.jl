# This file contains various methods to reduce systems of the PolSys form 
# The main idea is to have small functions that take a list of PolSys'es and 
# return a list with updated PolSys'es. 
# These can then be combined at will and fed to a function such as fixed_point 

export fixed_point

"""fixed_point(update_polsys::Function)::Function takes a function update_polsys that updates a list of PolSys'es and returns a 
function that keeps updating such a list until no changes have been found."""

function fixed_point(g::Function,eq=isequal)::Function 
    function h(x::Array{PolSys})::Array{PolSys}
        y = g(x)
        while !eq( y, x )
            x, y = y, g(x) 
        end
        return y
    end
    return h
end

export cleanup

"""cleanup(s::PolSys)::PolSys takes a PolSys and 
    (1) removes all duplicate polynomials, 
    (2) removes all nonzero variables from polynomials where they apear in each term,
    (3) removes all 0 polynomials, and 
    (4) removes all duplicate statements in the assumptions

   cleanup(ls::Array{PolSys})::Array{PolSys} maps cleanup to ls.
"""
function cleanup(systems::Array{PolSys})::Array{PolSys}
    map( cleanup, systems )
end

function cleanup( s::PolSys )::PolSys
end


export deduce_trivialities

"""deduce_trivialities( s::PolSys )::PolSys updates s by using binomials containing a constant term.
deduce_trivialities( systems::Array{PolSys} )::Array{Polsys} maps deduce_trivialities to systems.
"""
function deduce_trivialities( systems::Array{PolSys} )::Array{PolSys}
    map( deduce_trivialities, systems )
end

function deduce_trivialities( s::PolSys )::Array{PolSys}
end


export reduce_binomial_subsystem

"""reduce_binomial_subsystem(s::PolSys)::Array{PolSys} returns an array of PolSys' where 
    the binomial subsystem of s has been reduced such that the values of variables 
    that can be zero have been deduced and plugged in, and the subsystem of binomial equations
    have been upper triangularized. 

    Note: the list contains one reduced system per allowed configuration of variables that can be zero. 
   
    reduce_binomial_subsystem(ls::Array{PolSys})::Array{PolSys} maps reduce_binomials to ls.
"""

function reduce_binomial_subsystem(systems::Array{PolSys})::Array{PolSys}
    map( reduce_by_binomial_subsystem, systems )
end

function reduce_binomial_subsystem(s::PolSys)::Array{PolSys}

end


export solve_binomial_subsystem

"""reduce_binomials(s::PolSys)::Array{PolSys} returns an array of PolSys' where 
    the binomial subsystem of s has been reduced such that the values of variables 
    that can be zero have been deduced and plugged in, and the subsystem of binomial equations
    have been upper triangularized. 

    Note: the list contains one reduced system per allowed configuration of variables that can be zero. 
   
    reduce_binomials(ls::Array{PolSys})::Array{PolSys} maps reduce_binomials to ls.
"""

function solve_binomial_subsystem(systems::Array{PolSys})::Array{PolSys}
    map( reduce_by_binomials, systems )
end

function solve_binomial_subsystem(s::PolSys)::Array{PolSys}

end