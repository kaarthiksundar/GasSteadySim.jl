# function find_ub(ss::SteadySimulator, val::Float64, ub::Float64)::Float64
#     @assert ub > 0
#     while get_potential(ss, ub) < val
#         ub = 1.5 * ub
#     end 
#     return ub
# end 

# function find_lb(ss::SteadySimulator, val::Float64, lb::Float64)::Float64
#     @assert lb < 0
#     while get_potential(ss, lb) > val
#         lb = 1.5 * lb
#     end 
#     return lb
# end 

# function bisect(ss::SteadySimulator, lb::Float64, ub::Float64, val::Float64)::Float64  
#     @assert ub > lb
#     mb = 1.0
#     while (ub - lb) > TOL
#         mb = (ub + lb) / 2.0
#         if get_potential(ss, mb) > val
#             ub = mb
#         else
#             lb = mb 
#         end
#     end
#     return mb
# end

# invert_positive_potential(ss::SteadySimulator, val::Float64) = bisect(ss, 0.0, find_ub(ss, val, 1.0), val)

function calculate_slack_withdrawal(ss::SteadySimulator, id::Int, x_dof::Array)::Float64
    slack_withdrawal = 0.0
    for i in ref(ss, :incoming_dofs)[id]
        slack_withdrawal += x_dof[i]
    end 
    for i in ref(ss, :outgoing_dofs)[id]
        slack_withdrawal -= x_dof[i]
    end 
    return slack_withdrawal
end 


function update_solution_fields_in_ref!(ss::SteadySimulator, x_dof::Array)::NamedTuple
    flow_direction = true
    negative_flow_in_compressors = Int[]
    negative_nodal_pressures = Int[]

    for i in 1:length(x_dof)
        sym, local_id = ref(ss, :dof, i)
        if sym == :node
            
            if ref(ss, sym, local_id, "is_slack") == 1 
                ref(ss, sym, local_id)["withdrawal"] = calculate_slack_withdrawal(ss, local_id, x_dof)
            end 

            p_val = dof_to_pressure(x_dof[i]) 

            if (p_val < 0)
                push!(negative_nodal_pressures, local_id)
                ref(ss, sym, local_id)["pressure"] = p_val 
                ref(ss, sym, local_id)["density"] = get_density(ss, p_val)
            else  
                ref(ss, sym, local_id)["pressure"] = p_val
                ref(ss, sym, local_id)["density"] = get_density(ss, p_val)
            end
        end

        if sym == :pipe
            ref(ss, sym, local_id)["flow"] = x_dof[i]
            (x_dof[i] < 0) && (flow_direction = false)
        end

        if sym == :compressor
            ref(ss, sym, local_id)["flow"] = x_dof[i]
            if x_dof[i] < 0 && ref(ss, sym, local_id)["c_ratio"] > 1.0
                push!(negative_flow_in_compressors, local_id)
            end
        end

        (sym == :control_valve) && (ref(ss, sym, local_id)["flow"] = x_dof[i])
        (sym == :valve) && (ref(ss, sym, local_id)["flow"] = x_dof[i])
        (sym == :resistor) && (ref(ss, sym, local_id)["flow"] = x_dof[i])
        (sym == :loss_resistor) && (ref(ss, sym, local_id)["flow"] = x_dof[i])
        (sym == :short_pipe) && (ref(ss, sym, local_id)["flow"] = x_dof[i])
    end

    return (
            pipe_flow_dir = flow_direction, 
            compressors_with_neg_flow = negative_flow_in_compressors, 
            nodes_with_negative_pressures = negative_nodal_pressures
        )
end


function populate_solution!(ss::SteadySimulator)
    sol = ss.sol
    units = params(ss, :units)
    bc = ss.boundary_conditions

    function pressure_convertor(pu) 
        (units == 0) && (return pu * nominal_values(ss, :pressure)) 
        return pascal_to_psi(pu * nominal_values(ss, :pressure))
    end 

    function mass_flow_convertor(pu)
        kgps_to_mmscfd = get_kgps_to_mmscfd_conversion_factor(params(ss))
        (units == 0) && (return pu * nominal_values(ss, :mass_flow)) 
        return pu * nominal_values(ss, :mass_flow) * kgps_to_mmscfd
    end 

    function density_convertor(pu)
        return pu * nominal_values(ss, :density)
    end 
    
    for i in collect(keys(ref(ss, :node)))     
        sol["nodal_density"][i] = density_convertor(get_density(ss, ref(ss, :node, i, "pressure")))
        sol["nodal_pressure"][i] = pressure_convertor(ref(ss, :node, i, "pressure"))
    end

    for i in collect(keys(ref(ss, :pipe)))
        sol["pipe_flow"][i] = mass_flow_convertor(ref(ss, :pipe, i, "flow"))
    end
    
    if haskey(ref(ss), :compressor)
        for i in collect(keys(ref(ss, :compressor)))
            sol["compressor_flow"][i] = mass_flow_convertor(ref(ss, :compressor, i, "flow"))
        end
    end

    if haskey(ref(ss), :control_valve)
        for i in collect(keys(ref(ss, :control_valve)))
            sol["control_valve_flow"][i] = mass_flow_convertor(ref(ss, :control_valve, i, "flow"))
        end 
    end 

    if haskey(ref(ss), :control_valve)
        for i in bc[:control_valve_status][:off]
            sol["control_valve_flow"][i] = 0.0
        end 
    end

    if haskey(ref(ss), :valve)
        for i in collect(keys(ref(ss, :valve)))
            sol["valve_flow"][i] = mass_flow_convertor(ref(ss, :valve, i, "flow"))
        end 
    end 

    if haskey(ref(ss), :valve)
        for i in bc[:valve_status][:off]
            sol["valve_flow"][i] = 0.0
        end 
    end

    if haskey(ref(ss), :resistor)
        for i in collect(keys(ref(ss, :resistor)))
            sol["resistor_flow"][i] = mass_flow_convertor(ref(ss, :resistor, i, "flow"))
        end 
    end 

    if haskey(ref(ss), :loss_resistor)
        for i in collect(keys(ref(ss, :loss_resistor)))
            sol["resistor_flow"][i] = mass_flow_convertor(ref(ss, :loss_resistor, i, "flow"))
        end 
    end 

    if haskey(ref(ss), :short_pipe)
        for i in collect(keys(ref(ss, :short_pipe)))
            sol["short_pipe_flow"][i] = mass_flow_convertor(ref(ss, :short_pipe, i, "flow"))
        end 
    end 
    return
end 