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
                    for i ∈ 1:m1, j ∈ 1:m2, k ∈ 1:m3, l ∈ 1:m4
                        # [ F^{abc}_d ]^{(e,i,j)}_{(f,k,l)}
                        push!(
                            lbls,
                            [ a, b, c, d, e, i, j, f, k, l ]
                        )
                    end
                end
            end
        end
    end

    lbls
end

"""F_symbols( fr::FusionRing, field=ℚab) returns a tuple ( R, fsymb )
where R is a polynomial ring in the F-symbols fsymb over the requested field.
"""
function F_symbols( fr::FusionRing, field=ℚab)
    polynomial_ring( field, :ℱ => F_labels(fr), cached=true )
end



#TODO: figure out how to obtain indices of an F-symbol. Could be that its impossible since the indices might be part of a naming scheme and not properly stored :/
#
