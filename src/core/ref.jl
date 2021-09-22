function _add_components_to_ref!(ref::Dict{Symbol,Any}, data::Dict{String,Any})

    for (i, node) in get(data, "nodes", [])
        name = :node
        (!haskey(ref, name)) && (ref[name] = Dict())
        id = parse(Int64, i)
        ref[name][id] = Dict()
        @assert id == node["node_id"]
        ref[name][id]["id"] = id
        ref[name][id]["is_slack"] = node["slack_bool"]
        ref[name][id]["pressure"] = NaN
        ref[name][id]["density"] = NaN 
        ref[name][id]["withdrawal"] = NaN
        ref[name][id]["min_pressure"] = node["min_pressure"]
        ref[name][id]["max_pressure"] = node["max_pressure"]
    end

    for (i, pipe) in get(data, "pipes", [])
        name = :pipe
        (!haskey(ref, name)) && (ref[name] = Dict())
        id = parse(Int64, i)
        ref[name][id] = Dict()
        @assert id == pipe["pipe_id"]
        ref[name][id]["id"] = id
        ref[name][id]["fr_node"] = pipe["from_node"]
        ref[name][id]["to_node"] = pipe["to_node"]
        ref[name][id]["diameter"] = pipe["diameter"]
        ref[name][id]["area"] = pipe["area"]
        ref[name][id]["length"] = pipe["length"]
        ref[name][id]["friction_factor"] = pipe["friction_factor"]
        ref[name][id]["mass_flow"] = NaN
    end

    for (i, compressor) in get(data, "compressors", [])
        name = :compressor
        (!haskey(ref, name)) && (ref[name] = Dict())
        id = parse(Int64, i)
        ref[name][id] = Dict()
        @assert id == compressor["comp_id"]
        ref[name][id]["id"] = id
        ref[name][id]["to_node"] = compressor["to_node"]
        ref[name][id]["fr_node"] = compressor["from_node"]
        ref[name][id]["control_type"] = unknown_control
        ref[name][id]["c_ratio"] = NaN
        ref[name][id]["discharge_pressure"] = NaN
        ref[name][id]["flow"] = NaN
    end

    for (i, control_valve) in get(data, "control_valves", [])
        name = :control_valve
        (!haskey(ref, name)) && (ref[name] = Dict())
        id = parse(Int64, i)
        ref[name][id] = Dict()
        @assert id == control_valve["control_valve_id"]
        ref[name][id]["id"] = id
        ref[name][id]["to_node"] = control_valve["to_node"]
        ref[name][id]["fr_node"] = control_valve["from_node"]
        ref[name][id]["control_type"] = unknown_control
        ref[name][id]["c_ratio"] = NaN
        ref[name][id]["discharge_pressure"] = NaN
        ref[name][id]["flow"] = NaN
    end 
    return
end

function _add_index_info!(ref::Dict{Symbol, Any}, data::Dict{String, Any})
    dofid = 1
    ref[:dof] = Dict{Int64, Any}()
    
    for (i, node) in ref[:node]
        node[:dof] = dofid
        ref[:dof][dofid] = (:node, i)
        dofid += 1
    end

    for (i, pipe) in ref[:pipe]
        pipe[:dof] = dofid
        ref[:dof][dofid] = (:pipe, i)
        dofid += 1
    end

    for (i, compressor) in get(ref, :compressor, [])
        compressor[:dof] = dofid
        ref[:dof][dofid] = (:compressor, i)
        dofid += 1
    end

    for (i, control_valve) in get(ref, :control_valve, [])
        control_valve[:dof] = dofid
        ref[:dof][dofid] = (:control_valve, i)
        dofid += 1
    end
end

function _add_incident_dofs_info_at_nodes!(ref::Dict{Symbol,Any}, data::Dict{String,Any})
    ref[:incoming_dofs] = Dict{Int64, Vector{Int64}}()
    ref[:outgoing_dofs] = Dict{Int64, Vector{Int64}}()

    for (i, _) in ref[:node]
        ref[:incoming_dofs][i] = []
        ref[:outgoing_dofs][i] = []
    end

    for (_, pipe) in ref[:pipe]
        push!(ref[:incoming_dofs][pipe["to_node"]], pipe[:dof])
        push!(ref[:outgoing_dofs][pipe["fr_node"]], pipe[:dof])
    end

    for (_, compressor) in get(ref, :compressor, [])
        push!(ref[:incoming_dofs][compressor["to_node"]], compressor[:dof])
        push!(ref[:outgoing_dofs][compressor["fr_node"]], compressor[:dof])
    end

    for (_, control_valve) in get(ref, :control_valve, [])
        push!(ref[:incoming_dofs][control_valve["to_node"]], control_valve[:dof])
        push!(ref[:outgoing_dofs][control_valve["fr_node"]], control_valve[:dof])
    end

    return
end

function _add_pipe_info_at_nodes!(ref::Dict{Symbol,Any}, data::Dict{String,Any})
    ref[:incoming_pipes] = Dict{Int64, Vector{Int64}}()
    ref[:outgoing_pipes] = Dict{Int64, Vector{Int64}}()

    for (i, _) in ref[:node]
        ref[:incoming_pipes][i] = []
        ref[:outgoing_pipes][i] = []
    end

    for (id, pipe) in ref[:pipe]
        push!(ref[:incoming_pipes][pipe["to_node"]], id)
        push!(ref[:outgoing_pipes][pipe["fr_node"]], id)
    end
    return
end

function _add_compressor_info_at_nodes!(ref::Dict{Symbol,Any}, data::Dict{String,Any})
    ref[:incoming_compressors] = Dict{Int64, Vector{Int64}}()
    ref[:outgoing_compressors] = Dict{Int64, Vector{Int64}}()
    
    for (i, _) in ref[:node]
        ref[:incoming_compressors][i] = []
        ref[:outgoing_compressors][i] = []
    end

    for (id, compressor) in get(ref, :compressor, [])
        push!(ref[:incoming_compressors][compressor["to_node"]], id)
        push!(ref[:outgoing_compressors][compressor["fr_node"]], id)
    end
    return
end

function _add_control_valve_info_at_nodes!(ref::Dict{Symbol,Any}, data::Dict{String,Any})
    ref[:incoming_control_valves] = Dict{Int64, Vector{Int64}}()
    ref[:outgoing_control_valves] = Dict{Int64, Vector{Int64}}()
    
    for (i, _) in ref[:node]
        ref[:incoming_control_valves][i] = []
        ref[:outgoing_control_valves][i] = []
    end

    for (id, control_valve) in get(ref, :control_valve, [])
        push!(ref[:incoming_control_valves][control_valve["to_node"]], id)
        push!(ref[:outgoing_control_valves][control_valve["fr_node"]], id)
    end
    return
end


function build_ref(data::Dict{String,Any};
    ref_extensions=[])::Dict{Symbol,Any}

    ref = Dict{Symbol,Any}()
    _add_components_to_ref!(ref, data)

    for extension in ref_extensions
        extension(ref, data)
    end

    return ref
end