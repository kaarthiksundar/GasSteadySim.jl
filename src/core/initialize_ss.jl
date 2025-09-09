function initialize_simulator(data_folder::AbstractString;
    case_name::AbstractString="", 
    case_types::Vector{Symbol}=Symbol[],
    initial_guess_filename::AbstractString="",
    kwargs...)::SteadySimulator
    data = _parse_data(data_folder; 
        case_name=case_name, 
        case_types=case_types, 
        initial_guess_filename=initial_guess_filename
    )
    return initialize_simulator(data; kwargs...)
end

function initialize_simulator(data::Dict{String,Any}; 
    eos::Symbol=:ideal,
    use_potential_formulation::Bool=false, 
    potential_ratio_coefficients::Vector{Float64}=[0.0, 0.0, 1.0, 0.0])::SteadySimulator
    params, nominal_values = process_data!(data)
    make_per_unit!(data, params, nominal_values)
    bc = _build_bc(data)

    ref = build_ref(data, bc, ref_extensions= [
        _add_pipe_info_at_nodes!,
        _add_compressor_info_at_nodes!,
        _add_control_valve_info_at_nodes!,
        _add_valve_info_at_nodes!,
        _add_resistor_info_at_nodes!,
        _add_loss_resistor_info_at_nodes!,
        _add_short_pipe_info_at_nodes!,
        _add_index_info!,
        _add_incident_dofs_info_at_nodes!, 
        _add_pressure_node_flag!
        ]
    )
    
    # irrespective of the use_potential_formulation variable, ideal EoS always uses potentials 
    (eos == :ideal) && (_update_node_flag!(ref))
    # if eos is not ideal, then based on the flag call the node flag update function
    (use_potential_formulation == true) && (_update_node_flag!(ref))

    ig = _build_ig(data) 

    ss = SteadySimulator(data,
        ref,
        _initialize_solution(data),
        nominal_values,
        params,
        ig, bc,
        potential_ratio_coefficients,
        _get_eos(eos)...
    )

    # check if network is connected
    num_connected_components = check_network_connectivity(ss)
    @assert (num_connected_components == 1)  "Disconnected network ! $num_connected_components connected components found!"

    # _add_flow_bounds_to_ref!(ss)

    return ss
end

function check_network_connectivity(ss::SteadySimulator)::Int8
    V = ref(ss, :node) |> collect 
    new_node_from_old = Dict(V[i][1] => i for i in range(1, length(V)))
    old_node_from_new =[V[i][1] for i in range(1, length(V)) ]
    edge_component_universe = Vector{Symbol}([:pipe, :compressor, :control_valve, :valve, :short_pipe, :resistor, :loss_resistor])
    E = []
    for comp in edge_component_universe
        append!(E, get(ref(ss), comp, []))
    end
    g = SimpleGraph(length(V)) 

    for (i, edge) in E 
        fr = edge["fr_node"]
        to = edge["to_node"]
        if has_edge(g, new_node_from_old[fr], new_node_from_old[to])
            @warn "MORE THAN ONE EDGE BETWEEN NODES $(fr) AND $(to)"
            continue
        end
        add_edge!(g, new_node_from_old[fr], new_node_from_old[to])
    end 

    C_list = connected_components(g)
    
    return length(C_list)
end
