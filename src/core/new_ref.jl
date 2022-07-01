function _add_components_to_new_ref!(new_ref::Dict{Symbol,Any}, ref::Dict{Symbol,Any}, bc::Dict{Symbol,Any})

    new_ref[:pipe] = Dict()
    new_ref[:node] = Dict()

    for (key, pipe) in get(ref, :pipe, [])
        new_ref[:pipe][key] = Dict()
        new_ref[:pipe][key]["id"] = key
        new_ref[:pipe][key]["fr_node"] = pipe["fr_node"]
        new_ref[:pipe][key]["to_node"] = pipe["to_node"]
        new_ref[:pipe][key]["fr_factor"] = 1.0
        new_ref[:pipe][key]["to_factor"] = 1.0
        new_ref[:pipe][key]["diameter"] = pipe["diameter"]
        new_ref[:pipe][key]["area"] = pipe["area"]
        new_ref[:pipe][key]["length"] = pipe["length"]
        new_ref[:pipe][key]["friction_factor"] = pipe["friction_factor"]
        new_ref[:pipe][key]["flow"] = NaN
    end

    for (key, node) in get(ref, :node, [])
        new_ref[:node][key] = Dict()
        new_ref[:node][key]["id"] = key
        new_ref[:node][key]["is_slack"] = ref[:node][key]["is_slack"]
        new_ref[:node][key]["pressure"] = NaN
        new_ref[:node][key]["density"] = NaN 
        new_ref[:node][key]["withdrawal"] = 0.0
        new_ref[:node][key]["potential"] = NaN

        if haskey(bc[:node], key)
            if new_ref[:node][key]["is_slack"] == 1 #slack
                new_ref[:node][key]["withdrawal"] = NaN
                new_ref[:node][key]["pressure"] = bc[:node][key]["val"]
            else
                new_ref[:node][key]["withdrawal"] = bc[:node][key]["val"]
            end
        end

    end
    
    return
end

function _eliminate_compressor_edges_vertices!(new_ref::Dict{Symbol,Any}, ref::Dict{Symbol,Any}, bc::Dict{Symbol,Any})

    
    # (!haskey(ref(ss), :compressor)) && (return)
    for (comp_id, comp) in get(ref, :compressor, [])
        to_node = comp["to_node"]
        fr_node = comp["fr_node"]

        cmpr_val = bc[:compressor][comp_id]["val"]
        val1 = haskey(bc[:node], fr_node) ? bc[:node][fr_node]["val"] :   0.0
        val2 = haskey(bc[:node], to_node) ? bc[:node][to_node]["val"] : 0.0

        if new_ref[:node][to_node]["is_slack"] == 1
            new_ref[:node][fr_node]["is_slack"] = 1
            new_ref[:node][fr_node]["pressure"] = val2/cmpr_val
            println("Slack pressure moved from Node ", to_node, " to  Node", fr_node, "\n")
        elseif new_ref[:node][fr_node]["is_slack"] == 1
            new_ref[:node][fr_node]["pressure"] = val1
            println("Slack pressure at Node ", fr_node, "\n")

        else # both nodes slack makes system underdetermined
            new_ref[:node][fr_node]["withdrawal"] = val1 + val2
        end
        delete!(new_ref[:node], to_node)
    end

    new_ref[:incoming_pipes] =  Dict{Int64, Vector{Int64}}( i => [] for i in keys(new_ref[:node]) )
    new_ref[:outgoing_pipes] =  Dict{Int64, Vector{Int64}}( i => [] for i in keys(new_ref[:node]) )

    for (key, node) in get(new_ref, :node, [])
        new_ref[:incoming_pipes][key] = ref[:incoming_pipes][key]
        new_ref[:outgoing_pipes][key] = ref[:outgoing_pipes][key]
    end

    for (comp_id, comp) in get(ref, :compressor, [])
        edge_dof = comp["edge_dof"] 
        cmpr_val = bc[:compressor][comp_id]["val"]
        to_node = comp["to_node"]
        fr_node = comp["fr_node"]

        ip = ref[:incoming_pipes][to_node]
        for i in ip
            new_ref[:pipe][i]["to_node"] = fr_node
            new_ref[:pipe][i]["to_factor"] =  new_ref[:pipe][i]["to_factor"] * (cmpr_val^2)
            #add this pipe to list of incoming at vertex from_node
            push!(new_ref[:incoming_pipes][fr_node], i)
        end

        op = ref[:outgoing_pipes][to_node]
        for i in op
            new_ref[:pipe][i]["fr_node"] = fr_node
            new_ref[:pipe][i]["fr_factor"] =  new_ref[:pipe][i]["fr_factor"] * (cmpr_val^2)
            push!(new_ref[:outgoing_pipes][fr_node], i)
        end
    end

    return
end

function _add_index_info_new_ref!(new_ref::Dict{Symbol, Any}, ref::Dict{Symbol,Any}, bc::Dict{Symbol,Any})
    
    dofid = 0
    slack_dofid = 0
    vertex_dofid = 0
    edge_dofid = 0
    new_ref[:dof] = Dict{Int64, Any}()

    
    for (i, node) in new_ref[:node]
        
        if new_ref[:node][i]["is_slack"] == 1
            slack_dofid = slack_dofid - 1
            node["vertex_dof"] = slack_dofid
        else
            vertex_dofid += 1
            node["vertex_dof"] = vertex_dofid
            dofid += 1
            node["dof"] = dofid
            new_ref[:dof][dofid] = (:node, i)
        end
    end

    for (i, pipe) in new_ref[:pipe]
        dofid += 1
        edge_dofid += 1
        pipe["dof"] = dofid
        pipe["edge_dof"] = edge_dofid
        new_ref[:dof][dofid] = (:pipe, i)
    end


    new_ref[:total_slack_vertices] = abs(slack_dofid)
    new_ref[:total_nonslack_vertices] = vertex_dofid
    new_ref[:total_edges] = edge_dofid

    return

end




function _add_incident_dofs_info_at_nodes_new_ref!(new_ref::Dict{Symbol,Any}, ref::Dict{Symbol,Any}, bc::Dict{Symbol,Any})
    new_ref[:incoming_dofs] = Dict{Int64, Vector{Int64}}(
        i => [] for i in keys(new_ref[:node])
    )
    new_ref[:outgoing_dofs] = Dict{Int64, Vector{Int64}}(
        i => [] for i in keys(new_ref[:node])
    )

    for (_, pipe) in new_ref[:pipe]
        push!(new_ref[:incoming_dofs][pipe["to_node"]], pipe["dof"])
        push!(new_ref[:outgoing_dofs][pipe["fr_node"]], pipe["dof"])
    end

    # for (_, compressor) in get(ref, :compressor, [])
    #     push!(ref[:incoming_dofs][compressor["to_node"]], compressor["dof"])
    #     push!(ref[:outgoing_dofs][compressor["fr_node"]], compressor["dof"])
    # end

    

    return
end

# function _add_pipe_info_at_nodes_new_ref!(ref::Dict{Symbol,Any}, data::Dict{String,Any})
#     ref[:incoming_pipes] = Dict{Int64, Vector{Int64}}(
#         i => [] for i in keys(ref[:node])
#     )
#     ref[:outgoing_pipes] = Dict{Int64, Vector{Int64}}(
#         i => [] for i in keys(ref[:node])
#     )

#     for (id, pipe) in ref[:pipe]
#         push!(ref[:incoming_pipes][pipe["to_node"]], id)
#         push!(ref[:outgoing_pipes][pipe["fr_node"]], id)
#     end
#     return
# end

# function _add_compressor_info_at_nodes!(ref::Dict{Symbol,Any}, data::Dict{String,Any})
#     ref[:incoming_compressors] = Dict{Int64, Vector{Int64}}(
#         i => [] for i in keys(ref[:node])
#     )
#     ref[:outgoing_compressors] = Dict{Int64, Vector{Int64}}(
#         i => [] for i in keys(ref[:node])
#     )

#     for (id, compressor) in get(ref, :compressor, [])
#         push!(ref[:incoming_compressors][compressor["to_node"]], id)
#         push!(ref[:outgoing_compressors][compressor["fr_node"]], id)
#     end
#     return
# end


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


function build_new_ref(ref::Dict{Symbol,Any}, bc::Dict{Symbol,Any};
    ref_extensions=[])::Dict{Symbol,Any}

    new_ref = Dict{Symbol,Any}()
    _add_components_to_new_ref!(new_ref, ref, bc)

    for extension in ref_extensions
        extension(new_ref, ref, bc)
    end

    return new_ref
end