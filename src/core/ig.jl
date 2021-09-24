
function _build_ig(data::Dict{String,Any})::Dict{Symbol,Any}
    ig = Dict{Symbol,Any}()

    ig[:node] = Dict() 
    ig[:pipe] = Dict()
    ig[:compressor] = Dict()
    ig[:control_valve] = Dict()
    ig[:valve] = Dict()
    ig[:resistor] = Dict()
    ig[:loss_resistor] = Dict()
    ig[:short_pipe] = Dict()

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

    for (i, value) in get(data, "control_valve_flow", [])
        id = parse(Int64, i)
        ig[:control_valve][id] = value
    end 

    for (i, value) in get(data, "valve_flow", [])
        id = parse(Int64, i)
        ig[:valve][id] = value
    end 

    for (i, value) in get(data, "resistor_flow", [])
        id = parse(Int64, i)
        ig[:resistor][id] = value
    end 

    for (i, value) in get(data, "loss_resistor_flow", [])
        id = parse(Int64, i)
        ig[:loss_resistor][id] = value
    end 

    for (i, value) in get(data, "short_pipe_flow", [])
        id = parse(Int64, i)
        ig[:short_pipe][id] = value
    end 

    return ig
end 