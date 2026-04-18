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
    sort!(keys_in_order, by = x -> x[2])

    return [groups[key] for (key, _) in keys_in_order]
end