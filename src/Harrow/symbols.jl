# This file contains constructors for all kinds of symbols whose labels represent objects in a fusion category


export F_labels
# F_labels returns a list of 10-element lists corresponding to well-formed F-symbols
function F_labels(fr::FusionRing)::Vector{Vector{Int64}}
  r  = rank(fr)
  mt = multiplication_table(fr)
  sc = nzsc(fr)

  lbls = Vector{Int64}[]

  for lab1 in sc, lab2 in sc
    if lab1[3] == lab2[1]
      a, b, e = lab1;
      m1      = mt[a, b, e]
      c, d    = lab2[2:end];
      m2      = mt[e, c, d]
      for f in 1:r
        m3 = mt[b, c, f]
        m4 = mt[a, f, d]
        if m3 * m4 != 0
          for α in 1:m1, β in 1:m2, γ in 1:m3, δ in 1:m4
            push!(lbls, [a, b, c, d, e, α, β, f, γ, δ])
          end
        end
      end
    end
  end

  return lbls
end

export F_symbols

"""F_symbols( fr::FusionRing; ring=ℚ, symbol = :ℱ ) returns a couple ( R, dict )
where 
* R is a polynomial ring in the singly-indexed F-symbols [ symbol[1], ..., symbol[n] ] over the requested ring.
* dict maps the labels of the F-symbols to their singly-indexed representations
"""
function F_symbols(fr::FusionRing; ring = ℚ, symbol = :ℱ)
  fl = F_labels(fr)
  n  = length(fl)

  R, f = polynomial_ring(ring, symbol => 1:n; cached = true)
  dict = Dict(fl[i] => f[i] for i in 1:n)

  return (R, dict)
end

export vac_F_symbols

function vac_F_symbols( R::Ring, dict ) 
  return ( R, [ v for (k,v) in dict if 1 ∈ k[1:3] ] )
end

export F_matrices

function F_matrices( fr::FusionRing; ring = ℚ, symbol = :ℱ, triv_vac = true )
  R, dict = F_symbols( fr, ring = ring, symbol = symbol )
  F_matrices( R, dict, triv_vac = triv_vac )
end

# Takes polynomial ring in F-symbols and dict that translates
# labels of F-symbols to vars in ring and returns F-F_matrices
# in the F-symbols 
function F_matrices( R::Ring, dict; triv_vac = true )
  # group F-symbols by their first 4 labels
  labels = keys(dict)
  grouped_labels = gather_by( l -> l[1:4], labels )

  # function that maps labels to F-symbols
  function fs(v::Vector{Int64})
    if triv_vac
      return 1 ∈ v[1:3] ? R(1) : dict[v]
    else
      return  dict[v]
    end
  end

  # Create an F-matrix using a group of F-labels 
  function fmat( lbls::Vector{Vector{Int64}} )
    # associativity of fusion ring makes size of lbls square 
    m² = size(lbls,1) 
    m  = Int64(sqrt(m²))
    matrix( R, [ fs( lbls[ m*(i-1) + j ] ) for i in 1:m, j in 1:m  ] )
  end

  return ( R, dict, fmat.(grouped_labels) )

end
