
#= The idea of a tower is the following:
when solving a system of equations via backtracking it is generically interesting to 
start with the equations that contain few variables (or powerful equations, see later). 
Once a simplest equation has been identified, we can check it by looping over all values 
of its variables. When we choose the next simplest equation, we take account of the 
the fact that the values of variables for the previous equation(s) are known. 
Often, it happens that multiple equations are simultaneously checkable. We group these 
equations together. 

The weight function allows one to choose how to assess the complexity of a couple of
an equation and the unknown variables appearing in the equation.  

Example: if we want to build a tower based on the idea that the simplest equations are 
those with the least number of variables then we can represent the equations simply by 
couples (i, vars_i) where i is an integer and vars_i is the list of variables that appear 
in the i'th equation. The weight function would then be e -> size(e[2],1) and the 
tower would look like 
[ (vrs_1, ids_1 ), ..., (vrs_m, ids_m ) ]
where 
* ids_1 is the list of simplest equations assuming no values of the variables are known, and vrs_1 are the variables in those equations
* ids_2 is the list of the simplest equations assuming the values of vrs_1 are known, and vrs_2 are the remaining unknown variables in ids_2
* and so on

The function sort_and_group_by_knowns creates (what I call) a tower of couples of variables and lists of equations
=#


"""
sort_and_group_by_knowns( eqns, weight ) takes a list `eqns` of tuples ( eq_i, vrs_i ) that contain 
an equation (or identifier) eq_i and a list of variables apearing in that equation `vrs`. It returns a 
list of tuples ( vars_j, eqs_j ) where ids_j is a list of the equations that are 
fully determined via knowledge of vars_1, ..., vars_j. The order of the tuples is such that the weight of 
the elements of ids_j is the smallest out of all eq_i that are not contained in eqs_k (k<j) and for which 
the variables in vars_1, ... vars_{j-1} are not taken into account when calculating the weights.
"""
function sort_and_group_by_knowns( equations, weight::Function )
    newlist = []

    eqns = equations

    function removevars( vars, tuple )
        vrs, eqn = tuple 
        ( filter( x -> x ∉ vars, vrs ), eqn )
    end

    while !isempty(eqns)
        println(equations)
        # get index of eqn with lowest weight
        _, ind = findmin( weight.(eqns) ) 
        vars, eqn = eqns[ind]

        if isempty(vars) #  eqn doesn't contain new unknowns
            # append eqn to list of eqns of last tuple
            push!( newlist[end][2], eqn )
        else
            # add start new tuple
            push!( newlist, ( vars, [ eqn ] ) )
        end
       
        # remove known variables from other equations 
        eqns = map( eq -> removevars( vars, eq ), eqns )

        # remove equation
        deleteat!(eqns,ind)

    end

    newlist

end
