# abreviations for fields 
ℚ = QQ;
ℚb       = algebraic_closure(ℚ)
ℚab, ζ   = abelian_closure(ℚ; sparse = false)
sℚab, ζs = abelian_closure(ℚ; sparse = true)

function complement(a, b)
  ub      = unique(b)
  matches = []

  for x in a
    if x ∉ ub
      push!(matches, x)
    end
  end

  return matches
end

# List manipulation

function group_by(f, lis)
  isempty(lis) && return Dict()

  result = Dict{Any, Vector{eltype(lis)}}()
  for item in lis
    key = f(item)
    if !haskey(result, key)
      result[key] = eltype(lis)[]
    end
    push!(result[key], item)
  end
  return result
end

function gather_by(f, lis)
  isempty(lis) && return Vector{eltype(lis)}[]

  # Use vector of tuples to preserve order: 
  # (key, ind_of_first_occurrence)
  keys_in_order = Tuple{Any, Int}[]
  groups = Dict{Any, Vector{eltype(lis)}}()

  for (index, el) in enumerate(lis)
    key = f(el)
    # if key is new 
    if !haskey(groups, key)
      # create collection for key 
      groups[key] = eltype(lis)[]
      # store index of key to 
      push!(keys_in_order, (key, index))
    end
    # push element to group whose values evaluate to f(el)
    push!(groups[key], el)
  end

  # Sort keys_in_order by index
  sort!(keys_in_order; by = x -> x[2])

  return [groups[key] for (key, _) in keys_in_order]
end

# Sparse matrix functionality

# K is the ring to which elements of the vector belong. Necessary for
# when the vector is empty
function sparse_diagonal_matrix(K::Ring, v::Vector{T}) where {T}
  isempty(v) && return sparse_matrix(K)

  result = sparse_matrix(K)

  for i in 1:length(v)
    push!( result, sparse_row(K, [ ( i, v[i] ) ] ) )
  end

  return result
end




export to_composite_field

function to_composite_field(
  x::QQBarFieldElem; simplify_field = false, canonical_simplification = true
)
  arr, emb = to_composite_field([x]; simplify_field, canonical_simplification)
  return (arr[1], emb)
end

function to_composite_field(
  arr::Array{QQBarFieldElem}; simplify_field = false, canonical_simplification = true
)
  K, f = number_field(QQ, unique(arr))

  if simplify_field
    L, g = simplify(K; canonical = canonical_simplification)
    to_field_elem = x -> preimage(g, preimage(f, x))
    fg = hom(L, algebraic_closure(QQ), (f ∘ g ∘ gen)(L))
    return (to_field_elem.(arr), fg)
  else
    to_field_elem = x -> preimage(f, x)
    return (to_field_elem.(arr), f)
  end
end

export to_cyclotomic_field

function to_cyclotomic_field(
  x::QQBarFieldElem; simplify_field = false, canonical_simplification = true
)
  cfx, emb = to_composite_field(x; simplify_field, canonical_simplification)

  return to_cyclotomic_field(cfx, emb)
end

function to_cyclotomic_field(x::AbsSimpleNumFieldElem, emb)
  arr, emb, deg = to_cyclotomic_field([x], emb)

  return (arr[1], emb, deg)
end

function to_cyclotomic_field(arr::Array{AbsSimpleNumFieldElem}, emb)
  length(arr) === 0 && return (arr, emb)

  # Check parrent field of all fields are equal
  is_constant_array(parent.(arr)) || error("Elements of array should belong to same field")

  qqb = algebraic_closure(QQ)
  K   = parent(arr[1])
  C   = ray_class_field(K)
  deg = (Int ∘ minimum ∘ first ∘ conductor)(C)
  L,  = cyclotomic_field(deg)

  gen_K_as_cyclo = first(roots(L, defining_polynomial(K)))
  to_cyclo       = hom(K, L, gen_K_as_cyclo)

  rts = roots(qqb, defining_polynomial(L))
  for j in 1:deg
    emb_cyclo = hom(L, qqb, rts[j])

    if emb_cyclo(gen_K_as_cyclo) == emb(gen(K))
      return (to_cyclo.(arr), emb_cyclo, deg)
    else
      continue
    end
  end

  return error("Couldn't find embedding from cyclotomics into algebraic_closure(QQ)")
end

# combining two polynomial rings over the same base_field

# TODO: make itpossible to create common field from both rings if
# well-defined

"""merge(R1::Ring,R2::Ring)

takes two multivariate polynomial rings R1, R2, and returns a 4 tuple
R, newvars, i1, i2 where R is a polynomial ring in the combined number
of variables from R1, and R2, newvars is the list of new variables, and i1 resp
i2 are ring embeddings that map variables of R1 resp R2 to their counterparts
in R.
"""
function merge( R1::Ring, R2::Ring )
  base_ring(R1) != base_ring(R2) && error("Polynomial rings must have same base ring")

  s1, s2   = symbols(R1), symbols(R2)
  nv1, nv2 = number_of_variables(R1), number_of_variables(R2)

  R, newvars = polynomial_ring( base_ring(R1), vcat(s1,s2), cached = true  ) 

  i1 = hom( R1, R, newvars[1:nv1] )
  i2 = hom( R2, R, newvars[nv1+1:end] )

  R, newvars, i1, i2
end

"""merge(F1::AbsSimpleNumField, emb1, F2::AbsSimpleNumField, emb2; simplify = true )
returns a four-tuple K, ϕ, ι1, ι2 where
* K is the compositum field of the abstract fields F1 and F2
* ι1 resp ι2 are abstract embeddings of F1 resp F2 into K 
* ϕ is the embedding of K into QQBar such that ϕ ∘ ι1 = emb1 and ϕ ∘ ι2 = emb2,
  unless both emb1 == nothing and emb2 == nothing in which case ϕ will be nothing
if simplify = false, K will not be simplified
"""

function merge(F1::AbsSimpleNumField, emb1, F2::AbsSimpleNumField, emb2; simplify_field = true )
  if (emb1 === nothing && emb2 !== nothing) || (emb1 !== nothing && emb2 === nothing)
    error("Fields must either both be embedded or both be unembedded")
  end

  if F1 == QQabs && F2 == QQabs
    K   = QQabs
    ι1  = hom(F1,QQabs, QQabs(1))
    ι2  = hom(F2,QQabs, QQabs(1))
    emb = QQemb
    return K, emb, ι1, ι2
  end

  # TODO: implement other edge cases:
  # * F1 == F2 (with same and different embeddings),
  # * F1 ⊂ F2

  embedded = (emb1 !== nothing) && (emb2 !== nothing) 
    
  K, i1, i2 = compositum( F1, F2 )

  if embedded
    ϕ = hom(K, QQBar, emb1(gen(F1)) - emb2(gen(F2)) )

    if simplify_field
      S, is = simplify(K,canonical = true)

      ι1 = inv(is) ∘ i1 
      ι2 = inv(is) ∘ i2

      ϕs = ϕ ∘ is

      return S, ϕs, ι1, ι2
    else
      return K, ϕ, i1, i2
    end

  else
    #not embedded
    if simplify_field
      S, is = simplify(K,canonical = true)

      ι1 = inv(is) ∘ i1 
      ι2 = inv(is) ∘ i2

      return S, nothing, ι1, ι2
    else
      return K, nothing, i1, i2
    end
  end
end

#= tally counts the number of times an element apears in a vector and 
 returns a couple of vectors (els, counts) where els = unique 
 elements of v and counts = number of times they apear in v.
 If sort = false the elements in els are ordered by their first encounter 
 in the original original vector 
=#

function tally(v::AbstractVector; sort = false, sort_fun = identity)
  els = unique(v)

  counter(el) = count(==(el), v)

  counts = counter.(els)

  !sort && return (els, counts)

  s = sortperm(els; by = sort_fun)

  return (els[s], counts[s])
end

function is_constant_array(arr; equalfunc = ===)
  if isempty(arr)
    return true
  end
  first = arr[1]
  return all(equalfunc(element, first) for element in arr)
end
