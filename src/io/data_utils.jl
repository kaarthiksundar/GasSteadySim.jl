function _parse_data(data_folder::AbstractString; 
    case_name::AbstractString="", 
    case_types::Vector{Symbol}=Symbol[], 
    initial_guess_filename::AbstractString="")
    network_file = data_folder * "network.json"
    params_file = data_folder * "params"
    bc_file = data_folder * "bc"
    ig_file = data_folder * initial_guess_filename * ".json"
    if (:params in case_types)
        params_file = params_file * "_" * case_name * ".json"
    else 
        params_file = params_file * ".json"
    end 
    
    if (:bc in case_types)
        bc_file = bc_file * "_" * case_name * ".json"
    else 
        bc_file = bc_file * ".json"
    end 
    network_data = _parse_json(network_file)
    params_data = _parse_json(params_file)
    bc_data = _parse_json(bc_file)

    if isfile(ig_file)
        ig_data = _parse_json(ig_file)
        required_ig_fields = ["nodal_pressure", "pipe_flow", 
            "compressor_flow", "control_valve_flow", 
            "valve_flow", "resistor_flow", 
            "loss_resistor_flow", "short_pipe_flow"]
        filter!(p -> p.first in required_ig_fields, ig_data)
        data = merge(network_data, params_data, bc_data, ig_data)
    else
        @debug "initial guess file not provided, the code will generate an initial guess"
        data = merge(network_data, params_data, bc_data)
    end
    
    return data
end 

function _get_nominal_pressure(data::Dict{String,Any}, units)
    slack_pressures = []
    for (_, value) in get(data, "boundary_pslack", [])
        push!(slack_pressures, value)
    end 

    @assert length(slack_pressures) != 0
    (units == 1) && (return minimum(slack_pressures) *  6894.75729)
    return minimum(slack_pressures)
end 

function process_data!(data::Dict{String,Any})
    nominal_values = Dict{Symbol,Any}()
    params = Dict{Symbol,Any}()

    params_exhaustive = ["temperature", 
        "gas_specific_gravity",
        "specific_heat_capacity_ratio", 
        "nominal_length", 
        "nominal_velocity", 
        "nominal_pressure",
        "nominal_density",
        "units"]

    defaults_exhaustive = [288.706, 0.6, 1.4, 5000.0, NaN, NaN, NaN, 0]

    simulation_params = data["params"]
    
    key_map = Dict{String,String}()
    for k in keys(simulation_params)
        occursin("Temperature", k) && (key_map["temperature"] = k)
        occursin("Gas", k) && (key_map["gas_specific_gravity"] = k)
        occursin("Specific heat", k) &&
            (key_map["specific_heat_capacity_ratio"] = k)
        occursin("length", k) && (key_map["nominal_length"] = k)
        occursin("velocity", k) && (key_map["nominal_velocity"] = k)
        occursin("pressure", k) && (key_map["nominal_pressure"] = k)
        occursin("density", k) && (key_map["nominal_density"] = k)
        occursin("units", k) && (key_map["units"] = k)
    end

    # add "area" key to pipes in data
    for (_, pipe) in get(data, "pipes", [])
        pipe["area"] = pi * pipe["diameter"] * pipe["diameter"] * 0.25
    end

    # populating parameters
    for i in eachindex(params_exhaustive)
        param = params_exhaustive[i]
        default = defaults_exhaustive[i]
        if param == "units"
            if haskey(key_map, param)
                value = Int(simulation_params[key_map[param]])
                if (value == 0)
                    params[:units] = 0
                    params[:is_si_units] = 1
                    params[:is_english_units] = 0
                    params[:is_per_unit] = 0
                else
                    params[:units] = 1
                    params[:is_is_units] = 0
                    params[:is_english_units] = 1
                    params[:is_per_unit] = 0
                end
            else
                params[:units] = 0
                params[:is_si_units] = 1
                params[:is_english_units] = 0
                params[:is_per_unit] = 0
            end
            continue
        end
        
        key = get(key_map, param, false)
        if key != false
            value = simulation_params[key]
            params[Symbol(param)] = value
        else
            params[Symbol(param)] = default
        end
    end

    # other parameter calculations
    # universal gas constant (J/mol/K)
    params[:R] = 8.314
    # molecular mass of natural gas (kg/mol): M_g = M_a * G
    params[:gas_molar_mass] = 0.02896 * params[:gas_specific_gravity]
    params[:warning] = "R, temperature are in SI units. Rest are dimensionless"

    # sound speed (m/s): v = sqrt(R_g * T); 
    # R_g = R/M_g = R/M_a/G; R_g is specific gas constant; g-gas, a-air
    nominal_values[:sound_speed] = sqrt(params[:R] * params[:temperature] / params[:gas_molar_mass])
    nominal_values[:velocity] = params[:nominal_velocity] # choose based on mass flows
    nominal_values[:length] = params[:nominal_length]
    nominal_values[:area] = 1.0
    nominal_values[:pressure] = 
    if isnan(params[:nominal_pressure]) 
        _get_nominal_pressure(data, params[:units]) 
    else 
        params[:nominal_pressure]
    end
    nominal_values[:density] = 
    if isnan(params[:nominal_density]) 
        nominal_values[:pressure] / (nominal_values[:sound_speed]^2)
    else 
        params[:nominal_density]
    end 
    nominal_values[:velocity] =
    if isnan(params[:nominal_velocity])
       ceil(nominal_values[:sound_speed]/100)
    else
        params[:nominal_velocity]
    end
    nominal_values[:mass_flux] = nominal_values[:density] * nominal_values[:velocity]
    nominal_values[:mass_flow] = nominal_values[:mass_flux] * nominal_values[:area]
    nominal_values[:euler_num] = nominal_values[:pressure] / (nominal_values[:density] * nominal_values[:sound_speed]^2)
    nominal_values[:mach_num] = nominal_values[:velocity] / nominal_values[:sound_speed]
    
    return params, nominal_values
end