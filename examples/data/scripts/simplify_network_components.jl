
using GasSteadySim
using JSON

""" 
This function incorporates control_valves, valves, resistors, short_pipes etc. under the category of compressors
"""
function consolidate_all_non_pipe_edges_as_compressors(folder::AbstractString) 
    
    # folder = "./data/GasLib-2607/"
    network_file = folder * "network.json"
    bc_file = folder * "bc.json"
    params_file = folder * "params.json"

    
    new_network_filename = folder * "network_new.json"
    new_bc_filename =  folder * "bc_new.json"
    new_params_filename = folder * "params_new.json"

    params = GasSteadySim._parse_json(params_file)
    bc = GasSteadySim._parse_json(bc_file)
    network  = GasSteadySim._parse_json(network_file)
    num_compressors = length(network["compressors"])
    
    bc_new = Dict{String, Any}(
        "boundary_pslack" => Dict(), 
        "boundary_nonslack_flow" => Dict(), 
        "boundary_compressor" => Dict()
    )
    network_new = Dict{String, Any}(
        "nodes" => Dict(), 
        "pipes" => Dict(), 
        "compressors" => Dict()
    )

    network_new["nodes"] = network["nodes"]
    network_new["pipes"] = network["pipes"]
    network_new["compressors"] = network["compressors"]
    bc_new["boundary_pslack"] = bc["boundary_pslack"]
    bc_new["boundary_nonslack_flow"] = bc["boundary_nonslack_flow"]
    bc_new["boundary_compressor"] = bc["boundary_compressor"]

    compressor_like_components = Vector{String}(["control_valves", "valves"])

    new_id = num_compressors
    for (id_str, comp) in get(network, "control_valves", Dict())
        id_str in get(get(bc, "boundary_control_valve", Dict()), "off", []) && continue
        new_id += 1
        network_new["compressors"][string(new_id)] = network["control_valves"][id_str]
        network_new["compressors"][string(new_id)]["id"] = new_id #note this
        bc_new["boundary_compressor"][string(new_id)] = bc["boundary_control_valve"][id_str]
    end

    for (id_str, comp) in get(network, "valves", Dict())
        id_str in get(get(bc, "boundary_valve", Dict()), "off", []) && continue
        new_id += 1
        network_new["compressors"][string(new_id)] = network["valves"][id_str]
        network_new["compressors"][string(new_id)]["id"] = new_id
        bc_new["boundary_compressor"][string(new_id)] = Dict("control_type"=> 0,"value"=> 1.0)
    end

    pass_through_components = Vector{String}(["short_pipes", "resistors", "loss_resistors"])
    for item in pass_through_components
        for (id_str, comp) in get(network, item, Dict())
            new_id += 1
            network_new["compressors"][string(new_id)] = network[item][id_str]
            network_new["compressors"][string(new_id)]["id"] = new_id
            bc_new["boundary_compressor"][string(new_id)] = Dict("control_type"=> 0,"value"=> 1.0)
        end
    end

    open(new_params_filename, "w") do f 
        JSON.print(f, params, 2)
    end

    open(new_network_filename, "w") do f 
        JSON.print(f, network_new, 2)
    end

    open(new_bc_filename, "w") do f 
        JSON.print(f, bc_new, 2)
    end

    return
end 
