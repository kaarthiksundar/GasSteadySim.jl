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

function initialize_simulator(data::Dict{String,Any}; eos::Symbol=:ideal, 
    feasibility_model=:milp, num_partitions::Int=4)::SteadySimulator
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
        _add_incident_dofs_info_at_nodes!
        ]
    )

    ig = _build_ig(data) 

    ss = SteadySimulator(data,
        ref,
        _initialize_solution(data),
        nominal_values,
        params,
        ig, bc,
        JuMP.Model(),
        Dict(), Dict(),
        _get_eos(eos)...
    )

    construct_feasibility_model!(ss, 
        feasibility_model=feasibility_model, 
        num_partitions=num_partitions)
    return ss
end
