function _add_components_to_ref!(ref::Dict{Symbol,Any}, data::Dict{String,Any}, bc::Dict{Symbol,Any})

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
        ref[name][id]["potential"] = NaN
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
        ref[name][id]["flow"] = NaN
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

function _add_index_info!(ref::Dict{Symbol, Any}, data::Dict{String, Any})
    
    dofid = 0
    slack_dofid = 0
    vertex_dofid = 0
    edge_dofid = 0
    ref[:dof] = Dict{Int64, Any}()

    
    for (i, node) in ref[:node]
        

        if ref[:node][i]["is_slack"] == 1
            slack_dofid = slack_dofid - 1
            node["vertex_dof"] = slack_dofid
        else
            vertex_dofid += 1
            node["vertex_dof"] = vertex_dofid
            dofid += 1
            node["dof"] = dofid
            ref[:dof][dofid] = (:node, i)
        end
    end

    for (i, pipe) in ref[:pipe]
        dofid += 1
        edge_dofid += 1
        pipe["dof"] = dofid
        pipe["edge_dof"] = edge_dofid
        ref[:dof][dofid] = (:pipe, i)
    end

    for (i, compressor) in get(ref, :compressor, [])
        dofid += 1
        edge_dofid += 1
        compressor["dof"] = dofid
        compressor["edge_dof"] = edge_dofid
        ref[:dof][dofid] = (:compressor, i)
    end

    ref[:total_slack_vertices] = abs(slack_dofid)
    ref[:total_nonslack_vertices] = vertex_dofid
    ref[:total_edges] = edge_dofid

end

function _add_incident_dofs_info_at_nodes!(ref::Dict{Symbol,Any}, data::Dict{String,Any})
    ref[:incoming_dofs] = Dict{Int64, Vector{Int64}}(
        i => [] for i in keys(ref[:node])
    )
    ref[:outgoing_dofs] = Dict{Int64, Vector{Int64}}(
        i => [] for i in keys(ref[:node])
    )

    for (_, pipe) in ref[:pipe]
        push!(ref[:incoming_dofs][pipe["to_node"]], pipe["dof"])
        push!(ref[:outgoing_dofs][pipe["fr_node"]], pipe["dof"])
    end

    for (_, compressor) in get(ref, :compressor, [])
        push!(ref[:incoming_dofs][compressor["to_node"]], compressor["dof"])
        push!(ref[:outgoing_dofs][compressor["fr_node"]], compressor["dof"])
    end

    

    return
end

function _add_pipe_info_at_nodes!(ref::Dict{Symbol,Any}, data::Dict{String,Any})
    ref[:incoming_pipes] = Dict{Int64, Vector{Int64}}(
        i => [] for i in keys(ref[:node])
    )
    ref[:outgoing_pipes] = Dict{Int64, Vector{Int64}}(
        i => [] for i in keys(ref[:node])
    )

    for (id, pipe) in ref[:pipe]
        push!(ref[:incoming_pipes][pipe["to_node"]], id)
        push!(ref[:outgoing_pipes][pipe["fr_node"]], id)
    end
    return
end

function _add_compressor_info_at_nodes!(ref::Dict{Symbol,Any}, data::Dict{String,Any})
    ref[:incoming_compressors] = Dict{Int64, Vector{Int64}}(
        i => [] for i in keys(ref[:node])
    )
    ref[:outgoing_compressors] = Dict{Int64, Vector{Int64}}(
        i => [] for i in keys(ref[:node])
    )

    for (id, compressor) in get(ref, :compressor, [])
        push!(ref[:incoming_compressors][compressor["to_node"]], id)
        push!(ref[:outgoing_compressors][compressor["fr_node"]], id)
    end
    return
end


# function _add_pressure_node_flag!(ref::Dict{Symbol,Any}, data::Dict{String,Any})
#     ref[:is_pressure_node] = Dict{Int64, Bool}(
#         i => false for i in keys(ref[:node])
#     )

#     for (_, compressor) in get(ref, :compressor, [])
#         ref[:is_pressure_node][compressor["fr_node"]] = true 
#         ref[:is_pressure_node][compressor["to_node"]] = true
#     end 
# end 

# function _update_node_flag!(ref::Dict{Symbol,Any})
#     for i in keys(ref[:is_pressure_node])
#         ref[:is_pressure_node][i] = false
#     end 
# end 


function build_ref(data::Dict{String,Any}, bc::Dict{Symbol,Any};
    ref_extensions=[])::Dict{Symbol,Any}

    ref = Dict{Symbol,Any}()
    _add_components_to_ref!(ref, data, bc)

    for extension in ref_extensions
        extension(ref, data)
    end

    return ref
end