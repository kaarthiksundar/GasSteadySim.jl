
function _build_ig(data::Dict{String,Any})::Dict{Symbol,Any}
    ig = Dict{Symbol,Any}()

    ig[:node] = Dict() 
    ig[:pipe] = Dict()
    ig[:compressor] = Dict()
    ig[:regulator] = Dict()

    for (i, value) in get(data, "nodal_pressure", [])
        id = parse(Int64, i)
        ig[:node][id] = value
    end 

    for (i, value) in get(data, "pipe_flow", [])
        id = parse(Int64, i) 
        ig[:pipe][id] = value 
    end 

    for (i, value) in get(data, "compressor_flow", [])
        id = parse(Int64, i)
        ig[:compressor][id] = value
    end 

    for (i, value) in get(data, "regulator_flow", [])
        id = parse(Int64, i)
        ig[:regulator][id] = value
    end 

    return ig
end 