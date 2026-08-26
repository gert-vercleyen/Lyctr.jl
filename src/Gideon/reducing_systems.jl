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

function reduce_binomial_subsystem(s::PolSys)::Array{PolSys} end

export solve_binomial_subsystem

"""reduce_binomials(s::PolSys)::Array{PolSys} returns an array of PolSys' where 
    the binomial subsystem of s has been reduced such that the values of variables 
    that can be zero have been deduced and plugged in, and the subsystem of binomial equations
    have been upper triangularized. 

    Note: the list contains one reduced system per allowed configuration of variables that can be zero. 
   
    reduce_binomials(ls::Array{PolSys})::Array{PolSys} maps reduce_binomials to ls.
"""
# wz stands for without zeros: we assume none of the vars in the binomials can be zero!!!
function solve_binomial_subsystem_wz(
  systems::Array{PolSys}; expand_field = true
)::Array{PolSys}
  return map(s -> solve_binomial_system_wz(s; expand_field = expand_field), systems)
end


function solve_binomial_subsystem_wz(s::PolSys; expand_field = true, symbol = :z)
  !expand_field && error(
    "expand_field = false: Solving binomial subsystems over a fixed field is not yet implemented",
  )

  ϕ  = embedding(s)
  ϕ === nothing && error("Method only defined for embedded subfields of QQBar at the moment")

  isbinomial(p) = length(p) === 2

  # Filter binomials from system 
  bins = filter(isbinomial, polynomials(s))

  # If no binomials in system, return system
  isempty(bins) && return s

  # convert binomial equations to sparse matrix 
  sm, rhs = sparse_matrix_and_rhs(bins)

  nc = number_of_columns(sm)
  nr = number_of_rows(sm)

  # decompose matrix 
  ( L, D, R ) = diagonalize( sm )
  # get rank
  r = nnz(D)

  exprhs = powerdot( L, rhs )

  # TODO: need some logging functionality to convince user that system is indeed
  # (in)consistent

  # Check consistency
  !all(==(1),exprhs[r+1:end]) && pol_sys([],constraint(QQ,[false]),variables(s))

  # create r × r matrix inverse of the diagonal matrix d
  invdmat = sparse_diagonal_matrix( QQ, [ 1//D[i,i] for i in 1:r ] )

  # split R mat into pieces for discrete and continuous parts of solution
  Zmat = sub( R, 1:nc, 1:r )
  Cmat = sub( R, 1:nc, r+1:nc )

  cbr    = change_base_ring
  Zspace = mul_sparse(cbr(QQ,Zmat),invdmat)

  # DISCRETE PART OF SOLUTION
  #
  # set up all integer vectors that lead to inequivalent solutions
  denoms = denominator.( Matrix(Zspace) )
  lcms   = [ lcm( denoms[:,j] ) for j in 1:r ]
  ranges = map( i -> range(1,i), lcms )
  # TODO: following code feels very clunky, this shouldn't be so ugly...
  tuples = Iterators.product(ranges...)
  Zvecs  = reshape( [ collect(t) for t in tuples ], length(tuples), 1 )


  # set up cyclotomic field to express discrete solutions
  cd   = Int64(lcm([ D[i,i] for i in 1:r ]))
  cpol = cyclotomic_polynomial(cd)
  C, c = number_field( cpol )
  ic   = hom( C, QQBar, QQBar( ζ(cd) ) )

  # combine the old base field with the cyclotomic field
  #
  # we need to convert QQ to an abstract field in order to combine them
  K = base_field(s)
  if K === QQ 
    K = QQabs
    rhs = QQabs.(rhs)
  end
  ik = embedding(s)
  # S is the new base_field of the new pol system,
  # ι is an embedding from S into QQBar
  # iK embeds K into S and iL embeds L into S
  S, ι, iK, iL = merge( K, ik, C, ic )
  

  # set up exponentials of integer vectors. We might need to
  # expand the field over which the polsys is defined
  # TODO: what if we're using an abstract field? What about finite fields?
  
  
  cyclovec = [ c^(ZZ(cd//D[i,i])) for i in 1:r ]
  function expZvec( m )
    vec = cyclovec .^ m
    powerdot( Zmat, vec )
  end

  # CONSTANT PART OF SOLUTION 
  sL     = cbr( QQ, sub( L, 1:r, 1:nr ) )
  expmat = cbr( ZZ, mul_sparse( Zspace, sL ) )
  v0     = powerdot( expmat, rhs )

  # CONTINUOUS PART OF SOLUTION
  # TODO: R contains too many variables. A lot of the z[i] won't appear
  # in actual solution: should remove 0 columns of Cmat and only have
  # a variable per non-zero column
  R, z = polynomial_ring( S, :z => 1:(nc-r)  ) 
  FF   = fraction_field(R) # since we need to have negative powers of vars
  Cvec = powerdot( Cmat, FF.(z) )

  # COMBINED SOLUTION
  function sol( m )
    [ FF.(iK.(v0[i])) * FF.(iL.( expZvec(m)[i] )) * Cvec[i] for i in 1:nc ]
  end

  newsystems = PolSys[]

  for m in Zvecs
    solution   = sol(m)
    evsol(pol) = evaluate(pol,solution)
    newknowns  = Dict( k => evsol(v) for (k,v) in known_values(s) )
    newpols    = R.(numerator.(evsol.(polynomials(s))))

    #TODO:
    # * UPDATE CONSTRAINT!
    push!(newsystems, pol_sys(newpols,constraint(s),newknowns; emb =ι ) )
  end
  
  return newsystems
end

"""powerdot( v::SRow, w::SRow ) where v and w
   are two vectors with length n, returns w[1]^v[1] * ⋯ * w[n]^v[n]
   powerdot( M::SMat, w::SMat ) where w
   is a vector of length n, and M an m × n matrix returns a vector v of
   length n with coefficients
   v[i] =  w[1]^m[i,1] * ⋯ * w[n]^m[i,n]
"""
# TODO: need safety checks that vectors have right size
function powerdot( v::SRow, w::AbstractVector )
  ind = v.pos
  res = parent(w[1])(1)
  for i in ind
    res *= w[i]^v[i]
  end
  return res
end

function powerdot( v::AbstractVector, w::AbstractVector )
  res = parent(w[1])(1)
  for i in ind
    res *= w[i]^v[i]
  end
  return res
end

function powerdot( m::SMat, w::AbstractVector )::AbstractVector
  [ powerdot(m[i],w) for i in 1:nrows(m) ]
end

function powerdot( m::AbstractMatrix, w::AbstractVector )::AbstractVector
  [ powerdot(m[i,:],w) for i in 1:nrows ]
end

export sparse_matrix_and_rhs

function sparse_matrix_and_rhs(binomials::AbstractVector; bin_check = true)
  if is_empty(binomials)
    return (sparse_matrix(ZZ), Any[])
  end

  sparse_mat = sparse_matrix(ZZ)
  rhs        = []

  for binom in binomials
    srow, div = sparse_bin_row_rhs(binom; bin_check = bin_check)
    push!(sparse_mat, srow)
    push!(rhs, div)
  end

  return (sparse_mat, rhs)
end

# for binomial bin = α 𝐱^𝐔 + β 𝐱^𝐕 returns ( srow, rhs ) where srow is a 
# sparse row = 𝐮 - 𝐯 and rhs = -β/α
function sparse_bin_row_rhs(bin; bin_check = true)
  bin_check && length(bin) != 2 && error("argument "*string(bin)*" is not a binomial")
  𝐔 = exponent_vector(term(bin, 1), 1)
  𝐕 = exponent_vector(term(bin, 2), 1)

  s𝐔 = sparse_row(matrix(ZZ, [𝐔]))
  s𝐕 = sparse_row(matrix(ZZ, [𝐕]))
  
  α = coeff(bin, 1)
  β = coeff(bin, 2)

  return (s𝐔 - s𝐕, -β/α)
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
