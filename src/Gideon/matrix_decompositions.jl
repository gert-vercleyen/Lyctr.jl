#= This is code is based on the code from Jacob Strietelmeier's
  adaptation of the algorithms of Oscar.
  the code is meant to diagonalize/uppertriangularize a sparse ZZ 
  matrix and return the transformation matrices if desired
  input: matrix A, 
  output: V, D, U such that V * A * U = D 
  V, U can be omited by setting "left", resp "right" to false
=#

#┏━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━┓
#┃                                   diagonalize                                   ┃
#┗━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━┛

export diagonalize

function diagonalize(A::SMat, left::Bool=true, right::Bool=true)
  # n % 2 == 1 means A is not transposed
  n = 1 

  # sparse n x n identity mat 
  s1(n::Int64) = identity_matrix( SMat, ZZ, n )

  # hnf with our defaults 
  HNF( m ) = Hecke.hnf_kannan_bachem( m, truncate=false, full_hnf=true, auto=false )
  HNFT!( m, i ) = hnf_with_transform( m, truncate=false, full_hnf=true, auto=false )

  if left # want left transformation matrix
    C = s1( nrows(A) )
  end

  if right # want right transformation matrix
    D = s1( ncols(A) )
  end

  while(!is_diagonal(A))
    if n % 2 == 1 # A is not transposed
      if !left
        HNF!(A)
      else # TODO: this is almost the same code as for n % 1 == 0. Write function for this 
        I = s1( nrows(A) )

        HNFT!(A,I) 

      	R = sparse_matrix(ZZ, nrows(I), ncols(C))

        for (i, row) in enumerate(I.rows)
          rR = row * C
          R[i] = rR
        end
        C = R
      end

      A = transpose(A)
      n = n + 1

    else # A is transposed
      if !right
        HNF!(A)
      else
        I = s1( ncols(A)) 

        HNFT!(A,I)

        R = sparse_matrix(ZZ, nrows(I), ncols(D))

        for (i, row) in enumerate(I.rows)
          rR = row * D
          R[i] = rR
        end
        D = R
      end

      A = transpose(A)
      n = n + 1

    end
  end

  if n % 2 == 0
    A = transpose(A)
  end

  left && right && return ( C, A, transpose(D) )

  left && return ( C, A )

  right && return ( A, transpose(D) )

  return A
end    

#┏━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━┓
#┃                                        snf                                      ┃
#┗━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━┛

# TODO: implement sparse snf with transforms 

function hnf_with_transform(A::SMat{ZZRingElem}; truncate::Bool = false, full_hnf::Bool = true, auto::Bool = false, limit::Int=typemax(Int))
  B, T = hnf_kannan_bachem(A, Val(true); truncate, full_hnf, auto, limit)
  I = identity_matrix( SMat, ZZ, nrows(A))

  for (i, TT) in enumerate(T)
    apply_left!(I, TT)
  end
  return ( B, I )
end



