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

function initialize_simulator(data::Dict{String,Any}; eos::Symbol=:ideal)::SteadySimulator
    params, nominal_values = process_data!(data)
    make_per_unit!(data, params, nominal_values)
    bc = _build_bc(data)

    ref = build_ref(data, bc, ref_extensions= [
        _add_pipe_info_at_nodes!,
        _add_compressor_info_at_nodes!,
        _add_index_info!,
        _add_incident_dofs_info_at_nodes!
        ]
    )

    check_network_topology(ref)
    
    #can only be ideal so commenting out next line
    # (eos == :ideal) && (_update_node_flag!(ref)) 

    new_ref = build_new_ref(ref, bc, ref_extensions= [
        _eliminate_compressor_edges_vertices!,
        _add_index_info_new_ref!,
        _add_incident_dofs_info_at_nodes_new_ref!
        ]
    )
    
    ig = _build_ig(data) 

    ss = SteadySimulator(data,
        ref,
        new_ref,
        _initialize_solution(data),
        nominal_values,
        params,
        ig, bc,
        _get_eos(eos)...
    )


    return ss
end
