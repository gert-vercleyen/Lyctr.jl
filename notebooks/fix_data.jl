using JSON, UUIDs
datadir = "/home/gert/Projects/Lyctr.jl/src/data/"
#catdata = JSON.parsefile(datadir*"fusion_categories/fusion_categories_ug1.json")
props   = JSON.parsefile(datadir*"fusion_categories/cat_props.json")

dt = catdata["data"]

kys = collect(keys(dt))
ids  = Dict( k => string(UUIDs.uuid1()) for k in kys )
ringid = FusionRings.uuid
fawc   = FusionRings.fawc

_nonames = Dict(
  "quantum_group_like" => missing,
  "group_like"         => missing,
  "physics"            => missing,
  "miscellaneous"      => missing,
)

function mscnames(v::Vector)
  return Dict(
    "quantum_group_like" => missing,
    "group_like"         => missing,
    "physics"            => missing,
    "miscellaneous"      => string.(v),
  )
end

function fix_data( k )
  d = Dict( ky => vl for (ky,vl) in dt[k] ) 

  code = parse.(Int64,split(k[6:end],"_"))
  scs  = k[6:end]

  fr       = fawc(code[1:4]...)
  ringuuid = ringid(fr)
  catuuid  = uuids[k]

  # set tex names
  tn = FusionRings.tex_names(fr)["miscellaneous"] 
  endcode = code[5:7]
  function tocatname(str::String)::String
    string("[",str,"]_{",join(string.(endcode),","),"}")
  end
  
  d["anyonwiki_code"] = code
  d["texnames"]    = mscnames(tocatname.(tn))
  d["is_modular"]  = props[scs]["modular"]
  d["is_ribbon"]   = props[scs]["ribbon"]
  d["is_spherical"]   = props[scs]["spherical"]
  d["fusion_ring"] = ringuuid
  d["uuid"]        = catuuid
  delete!(d,"id")
  d
end

fixed_data = [  fix_data(k) for k in kys ]
