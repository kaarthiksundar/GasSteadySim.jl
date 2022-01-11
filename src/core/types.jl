struct SteadySimulator
    data::Dict{String,Any}
    ref::Dict{Symbol,Any}
    sol::Dict{String,Any}
    nominal_values::Dict{Symbol,Any}
    params::Dict{Symbol,Any}
    initial_guess::Dict{Symbol,Any}
    boundary_conditions::Dict{Symbol,Any}
    pu_eos_coeffs::Function
    pu_pressure_to_pu_density::Function
    pu_density_to_pu_pressure::Function
end

ref(ss::SteadySimulator) = ss.ref
ref(ss::SteadySimulator, key::Symbol) = ss.ref[key]
ref(ss::SteadySimulator, key::Symbol, id::Int64) = ss.ref[key][id]
ref(ss::SteadySimulator, key::Symbol, id::Int64, field) = ss.ref[key][id][field]

params(ss::SteadySimulator) = ss.params
params(ss::SteadySimulator, key::Symbol) = ss.params[key]

nominal_values(ss::SteadySimulator) = ss.nominal_values
nominal_values(ss::SteadySimulator, key::Symbol) = ss.nominal_values[key]

initial_pipe_mass_flow(ss::SteadySimulator, id::Int64) = 
    ss.initial_guess[:pipe][id]

initial_compressor_flow(ss::SteadySimulator, id::Int64) = 
    ss.initial_guess[:compressor][id]

initial_nodal_pressure(ss::SteadySimulator, id::Int64) = 
    ss.initial_guess[:node][id]

initial_control_valve_flow(ss::SteadySimulator, id::Int64) = 
    ss.initial_guess[:control_valve][id]

initial_valve_flow(ss::SteadySimulator, id::Int64) = 
    ss.initial_guess[:valve][id]

initial_resistor_flow(ss::SteadySimulator, id::Int64) = 
    ss.initial_guess[:resistor][id]

initial_loss_resistor_flow(ss::SteadySimulator, id::Int64) = 
    ss.initial_guess[:loss_resistor][id]

initial_short_pipe_flow(ss::SteadySimulator, id::Int64) = 
    ss.initial_guess[:short_pipe][id]

function control(ss::SteadySimulator,
    key::Symbol, id::Int64)::Tuple{CONTROL_TYPE,Float64}
    (key == :node) && (return get_nodal_control(ss, id))
    (key == :compressor) && (return get_compressor_control(ss, id))
    (key == :control_valve) && (return get_control_valve_control(ss, id))
    @error "control available only for nodes, compressors, and control_valves"
    return CONTROL_TYPE::unknown_control, 0.0
end

get_eos_coeffs(ss::SteadySimulator) = ss.pu_eos_coeffs(nominal_values(ss), params(ss))
get_pressure(ss::SteadySimulator, density) = ss.pu_density_to_pu_pressure(density, nominal_values(ss), params(ss))
get_density(ss::SteadySimulator, pressure) = ss.pu_pressure_to_pu_density(pressure, nominal_values(ss), params(ss))

TOL = 1.0e-5

function get_nodal_control(ss::SteadySimulator,
    id::Int64)::Tuple{CONTROL_TYPE,Float64}
    if !haskey(ss.boundary_conditions[:node], id)
        return flow_control, 0.0
    end
    control_type = ss.boundary_conditions[:node][id]["control_type"]
    val = ss.boundary_conditions[:node][id]["val"]
    return control_type, val
end

function get_compressor_control(ss::SteadySimulator,
    id::Int64)::Tuple{CONTROL_TYPE,Float64}
    control_type = ss.boundary_conditions[:compressor][id]["control_type"]
    val = ss.boundary_conditions[:compressor][id]["val"]
    return CONTROL_TYPE(control_type), val
end

function get_control_valve_control(ss::SteadySimulator,
    id::Int64)::Tuple{CONTROL_TYPE,Float64}
    control_type = ss.boundary_conditions[:control_valve][id]["control_type"]
    val = ss.boundary_conditions[:control_valve][id]["val"]
    return CONTROL_TYPE(control_type), val
end

@enum CONTROL_TYPE begin
    c_ratio_control = 0
    discharge_pressure_control = 1
    flow_control = 2
    pressure_control = 3
    unknown_control = 100
end

@enum SOLVER_STATUS begin 
    successfull = 0 
    successfull_after_pressure_correction = 1
    initial_nl_solve_failure = 2
    pressure_correction_nl_solve_failure = 3
    pressure_correction_failure = 4
    compressor_flow_negative = 5
end

struct SolverReturn 
    status::SOLVER_STATUS
    iterations::Int 
    residual_norm::Float64 
    time::Float64 
    solution::Vector{Float64}
    negative_flow_in_compressors::Vector{Int64}
end 
