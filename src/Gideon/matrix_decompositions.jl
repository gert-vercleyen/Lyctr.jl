#= This is code is based on the code from Jacob Strietelmeier's
  adaptation of the algorithms of Oscar.
  the code is meant to diagonalize/uppertriangularize a sparse ZZ 
  matrix and return the transformation matrices if desired
  input: matrix A, 
  output: L, D, R such that L * A * R = D 
  L, R can be omitted (returning nothing instead) by setting "left", resp "right" to false
=#

#┏━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━┓
#┃                                   diagonalize                                   ┃
#┗━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━┛

export diagonalize

function diagonalize(A::SMat, left::Bool = true, right::Bool = true)
  # n % 2 == 1 means A is not transposed
  n = 1 

  # sparse n x n identity mat 
  s1(n::Int64) = identity_matrix(SMat, ZZ, n)

  # hnf with specific defaults 
  HNF(m) = Hecke.hnf_kannan_bachem(m, truncate = false, full_hnf = true, auto = false)
  HNFT(m, i) = hnf_with_transform(m, truncate = false, full_hnf = true, auto = false)

  # if want left transformation matrix
    C = left ? s1(nrows(A)) : nothing

  # if want right transformation matrix
    D = right ? s1(ncols(A)) : nothing

  nr = nrows(A)
  nc = ncols(A)

  while !is_diagonal(A)
    if n % 2 == 1 # A is not transposed
      if !left
        HNF!(A)
      else 
        I = s1(nr)

        A, I = HNFT(A, I)

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
        I = s1(nc) 

        A, I = HNFT(A, I)

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

  return (C, A, transpose(D))
end    

#┏━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━┓
#┃                                        snf                                      ┃
#┗━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━┛

function hnf_with_transform(
  A::SMat{ZZRingElem};
  truncate::Bool = false,
  full_hnf::Bool = true,
  auto::Bool = false,
  limit::Int = typemax(Int),
)
  B, T = Hecke.hnf_kannan_bachem(A, Val(true); truncate, full_hnf, auto, limit)

  I = identity_matrix(SMat, ZZ, nrows(A))

  for m in T
        Hecke.apply_left!(I,m)
  end

  return (B, I )
end
