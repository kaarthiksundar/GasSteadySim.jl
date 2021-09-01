function add_components_to_ref!(ref::Dict{Symbol,Any}, data::Dict{String,Any})

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
    return
end

function add_index_info!(ref::Dict{Symbol, Any}, data::Dict{String, Any})
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

    for (i, compressor) in ref[:compressor]
        compressor[:dof] = dofid
        ref[:dof][dofid] = (:compressor, i)
        dofid += 1
    end
end

function add_incident_dofs_info_at_nodes!(ref::Dict{Symbol,Any}, data::Dict{String,Any})

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

    return
end


function build_ref(data::Dict{String,Any};
    ref_extensions=[])::Dict{Symbol,Any}

    ref = Dict{Symbol,Any}()
    add_components_to_ref!(ref, data)

    for extension in ref_extensions
        extension(ref, data)
    end

    return ref
end