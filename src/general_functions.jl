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
