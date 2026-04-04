# This file contains constructors for all kinds of symbols whose labels represent objects in a fusion category

# Some abreviations for fields
ℚ        = QQ;
ℚb       = algebraic_closure( ℚ )
ℚab, ζ   = abelian_closure( ℚ, sparse=false )
sℚab, ζs = abelian_closure( ℚ, sparse=true )


# F_labels returns a list of 10-element lists corresponding to well-formed F-symbols
function F_labels(fr::FusionRing)
    r  = rank( fr )
    mt = multiplication_table( fr )
    sc = nzsc( fr )

    lbls = Vector{Int64}[]

    for lab1 in sc, lab2 in sc
        if lab1[3] == lab2[1]
            a, b, e = lab1;        m1 = mt[ a, b, e ]
            c, d    = lab2[2:end]; m2 = mt[ e, c, d ]
            for f in 1:r
                m3 = mt[ b, c, f ]
                m4 = mt[ a, f, d ]
                if m3 * m4 != 0
                    for α ∈ 1:m1, β ∈ 1:m2, γ ∈ 1:m3, δ ∈ 1:m4
                        # [ F^{abc}_d ]^{(e,i,j)}_{(f,k,l)}
                        push!(
                            lbls,
                            [ a, b, c, d, e, α, β, f, γ, δ ]
                        )
                    end
                end
            end
        end
    end

    lbls
end

"""F_symbols( fr::FusionRing, field=ℚab, symbol = :ℱ ) returns a fourtuple ( R, fvar, dict, idict )
where 
* R is a polynomial ring in the singly-indexed F-symbols fvar over the requested field.
* fvar is an array of singly-indexed F-symbols [ ℱ[1], ..., ℱ[n] ]
* dict maps the singly-index F-symbols to their labels (arrays with 10 integers)
* idict maps the labels of the F-symbols to their singly-indexed representations
"""
function F_symbols( fr::FusionRing; ring = ℚb, symbol = :ℱ )
    fl = F_labels(fr)
    n  = length(F_labels)

    R, f  = polynomial_ring( field, :ℱ => 1:n, cached=true )
    dict  = Dict( f[i] => fl[i] for i in 1:n )
    idict = Dict( fl[i] => f[i] for i in 1:n ) 

    ( R, f, dict, idict )
end



