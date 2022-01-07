function _build_bc(data::Dict{String,Any})::Dict{Symbol,Any}
    bc = Dict{Symbol,Any}()

    bc[:node] = Dict() 
    bc[:compressor] = Dict()
    bc[:control_valve] = Dict()
    bc[:control_valve_status] = Dict()
    bc[:valve_status] = Dict()

    if !haskey(data, "boundary_pslack") 
        @error "Simulator requires at least one slack node pressure for uniqueness of solutions"
    end 

    for (i, value) in get(data, "boundary_pslack", [])
        id = parse(Int64, i)
        bc[:node][id] = Dict( 
            "val" => value, 
            "control_type" => pressure_control
        )
    end 

    for (i, value) in get(data, "boundary_nonslack_flow", [])
        id = parse(Int64, i)
        bc[:node][id] = Dict(
            "val" => value,
            "control_type" => flow_control
        )
    end 

    for (i, value) in get(data, "boundary_compressor", [])
        id = parse(Int64, i)
        if value["control_type"] != 0 
            @error "Simulator does not support compressor flow or discharge pressure specification for compressor boundary condition"
        end 
        bc[:compressor][id] = Dict(
            "val" => value["value"],
            "control_type" => value["control_type"]
        )
    end 

    for (i, value) in get(data, "boundary_control_valve", [])
        (i == "on") && (bc[:control_valve_status][:on] = value; continue) 
        (i == "off") && (bc[:control_valve_status][:off] = value; continue) 
        id = parse(Int64, i)
        if value["control_type"] != 0 
            @error "Simulator does not support compressor flow or discharge pressure specification for control valve boundary condition"
        end 
        bc[:control_valve][id] = Dict(
            "val" => value["value"],
            "control_type" => value["control_type"]
        )
    end 

    for (i, value) in get(data, "boundary_valve", [])
        (i == "on") && (bc[:valve_status][:on] = value)
        (i == "off") && (bc[:valve_status][:off] = value)
    end

    return bc
end 