# This file contains various methods to reduce systems of the PolSys form 
# The main idea is to have small functions that take a list of PolSys'es and 
# return a list with updated PolSys'es. 
# These can then be combined at will and fed to a function such as fixed_point 

export fixed_point

"""fixed_point(update_polsys::Function)::Function takes a function update_polsys that updates a list of PolSys'es and returns a 
function that keeps updating such a list until no changes have been found."""

function fixed_point(g::Function, eq = isequal)::Function
  function h(x::Array{PolSys})::Array{PolSys}
    y = g(x)
    while !eq(y, x)
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
  return map(cleanup, systems)
end

function cleanup(s::PolSys)::PolSys end

export deduce_trivialities

"""deduce_trivialities( s::PolSys )::PolSys updates s by using binomials containing a constant term.
deduce_trivialities( systems::Array{PolSys} )::Array{Polsys} maps deduce_trivialities to systems.
"""
function deduce_trivialities(systems::Array{PolSys})::Array{PolSys}
  return map(deduce_trivialities, systems)
end

function deduce_trivialities(s::PolSys)::Array{PolSys} end

export reduce_binomial_subsystem

"""reduce_binomial_subsystem(s::PolSys)::Array{PolSys} returns an array of PolSys' where 
    the binomial subsystem of s has been reduced such that the values of variables 
    that can be zero have been deduced and plugged in, and the subsystem of binomial equations
    have been upper triangularized. 

    Note: the list contains one reduced system per allowed configuration of variables that can be zero. 
   
    reduce_binomial_subsystem(ls::Array{PolSys})::Array{PolSys} maps reduce_binomials to ls.
"""

function reduce_binomial_subsystem(systems::Array{PolSys})::Array{PolSys}
  return map(reduce_by_binomial_subsystem, systems)
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
# wz stands for without zeros: we assume none of the vars in the binomials can be zero!!!
function solve_binomial_system_wz( systems::Array{PolSys}; expand_field = true )::Array{PolSys}
  return map( s -> solve_binomial_system_wz( s, expand_field = expand_field ), systems)
end


function solve_binomial_system_wz(s::PolSys, expand_field = true )::Array{PolSys}
  !expand_field && error("expand_field = false: Solving binomial subsystems over a fixed field is not yet implemented")

  # Filter binomials from system 
  bins = filter( is_binomial , polynomials(s) )

  # convert binomial equations to sparse matrix 
  sm, rhs = sparse_matrix_and_rhs(bins)

  # decompose matrix 

  # get rank 

  # check if rank is 0

  #  

end

export sparse_matrix_and_rhs

function sparse_matrix_and_rhs( binomials::AbstractVector; bin_check = true ) 

  if is_empty(binomials)
    return ( sparse_matrix(ZZ), Any[] )
  end

  sparse_mat = sparse_matrix(ZZ)
  rhs        = []

  for binom in binomials
    srow, div = sparse_bin_row_rhs(binom, bin_check = bin_check )
    push!( sparse_mat, srow )
    push!( rhs, div )
  end

  ( sparse_mat, rhs )

end

# for binomial bin = α 𝐱^𝐔 + β 𝐱^𝐕 returns ( srow, rhs ) where srow is a 
# sparse row = 𝐮 - 𝐯 and rhs = -β/α
function sparse_bin_row_rhs( bin; bin_check = true )
  bin_check && !is_binomial(bin) && error("argument is not a binomial")
  𝐔 = exponent_vector( term( bin, 1 ), 1 )
  𝐕 = exponent_vector( term( bin, 2 ), 1 )

  s𝐔 = sparse_row( matrix( ZZ, [ 𝐔 ] ) )
  s𝐕 = sparse_row( matrix( ZZ, [ 𝐕 ] ) )
  
  α = coeff( bin, 1 )
  β = coeff( bin, 2 )

  ( s𝐔 - s𝐕, -β/α )
end

#= function most likely not needed unless converting to sparse 
   rows via matrix(ZZ,...) in sparse_bin_row_rhs is too costly
# return the index of the variables and their exponents
# for a term 
function pos_exponent(term)

  exponents = filter( !=(0) , exponent_vector(term,1))
  varind    = var_indices(term)

  zip( varind, exponents )
end
=#

#function reduce_semi_lin_mod_ℤ( sm::, rhs::, :z )
#end


