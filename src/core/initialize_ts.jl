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
    ref = build_ref(data, ref_extensions= [
        add_pipe_info_at_nodes!,
        add_compressor_info_at_nodes!,
        add_index_info!,
        add_incident_edge_dofs_info_at_nodes!
        ]
    )

    bc = build_bc(data)
    feasibility_model = JuMP.Model() # get_feasibility_model(ref, bc, eos)

    ss = SteadySimulator(data,
        ref,
        initialize_solution(data),
        nominal_values,
        params,
        build_ig(data),
        bc, 
        feasibility_model,
        get_eos(eos)...
    )
    return ss
end
