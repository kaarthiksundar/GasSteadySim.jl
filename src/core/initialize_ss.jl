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
    initial_guess_from_opt::Bool=true,
    optimizer=_highs_optimizer
    )::SteadySimulator
    params, nominal_values = process_data!(data)
    params[:initial_guess_from_opt] = initial_guess_from_opt
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
    
    (eos == :ideal) && (_update_node_flag!(ref))
    
    ig = _build_ig(data) 

    ss = SteadySimulator(data,
        ref,
        _initialize_solution(data),
        nominal_values,
        params,
        ig, bc, OptModel(optimizer),
        _get_eos(eos)...
    )

    _add_flow_bounds_to_ref!(ss)
    _populate_lp_model!(ss)

    if initial_guess_from_opt
        _solve_lp_model!(ss)
        if _is_optimal(ss)
            _populate_lp_solution!(ss)
        end 
    end 

    return ss
end
