
#┏━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━┓
#┃                   Properties of SkeletalFusionCats                              ┃
#┗━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━┛


# Fusion ring
export fusion_ring
fusion_ring(sfc::SkeletalFusionCat)::FusionRing = sfc.fr

# ID
export fusion_cat_id
fusion_cat_id(sfc::SkeletalFusionCat)::String = sfc.id

# F-symbols
export F_symbols
F_symbols(sfc::SkeletalFusionCat)::FSymbols = sfc.fsymb

# R-symbols
export R_symbols
R_symbols(sfc::SkeletalFusionCat)::Union{RSymbols,Missing,Nothing} = sfc.rsymb

# P-symbols
export P_symbols
P_symbols(sfc::SkeletalFusionCat)::Union{PSymbols,Missing} = sfc.psymb

# Base field
export base_field
base_field(sfc::SkeletalFusionCat) = sfc.basefield

# Embedding
export field_embedding
field_embedding(sfc::SkeletalFusionCat) = sfc.embedding

# In minimal field
export is_in_min_field
is_in_min_field(sfc::SkeletalFusionCat)::Bool = sfc.inminfield

# Minimal fields
export minimal_fields
minimal_fields(sfc::SkeletalFusionCat) = sfc.minimalfields

# In unitary gauge
export is_in_unitary_gauge
is_in_unitary_gauge(sfc::SkeletalFusionCat)::Bool = sfc.inunitarygauge

# Names
export names
names(sfc::SkeletalFusionCat)::Vector{String} = sfc.nms

# TeX names
export tex_names
tex_names(sfc::SkeletalFusionCat)::Vector{String} = sfc.texnms

# Inverse gauge-split transform
export inverse_gauge_split_transform
inverse_gauge_split_transform(sfc::SkeletalFusionCat) = sfc.igst

# Gauge-split basis
export gauge_split_basis
gauge_split_basis(sfc::SkeletalFusionCat) = sfc.gsb

# Realizations
export realizations
realizations(sfc::SkeletalFusionCat) = sfc.rlztns

# Is pivotal
export is_pivotal
is_pivotal(sfc::SkeletalFusionCat)::Union{Missing,Bool} = sfc.isp

# Is unitary
export is_unitary
is_unitary(sfc::SkeletalFusionCat)::Union{Missing,Bool} = sfc.isu

# Is braided
export is_braided
is_braided(sfc::SkeletalFusionCat)::Union{Missing,Bool} = sfc.isb

# Is modular
export is_modular
is_modular(sfc::SkeletalFusionCat)::Union{Missing,Bool} = sfc.ism



# polynomial ring in all symbols including gauge symbols
function polynomial_ring(cat::SkeletalFusionCat)::Ring end
