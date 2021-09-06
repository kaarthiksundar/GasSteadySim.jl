
function update_solution_fields_in_ref!(ss::SteadySimulator, x_dof::Array)

    flow_direction = true
    compressor_flow_direction = true
    for i = 1: length(x_dof)
        sym, local_id = ref(ss, :dof, i)
        if sym == :node
            ref(ss, sym, local_id)["pressure"] = x_dof[i]
            ref(ss, sym, local_id)["density"] = get_density(ss, x_dof[i])
            ctrl_type, val = control(ss, :node, local_id)
            if ctrl_type == flow_control
                ref(ss, sym, local_id)["withdrawal"] = val
            end
        end

        if sym == :pipe
            ref(ss, sym, local_id)["mass_flow"] = x_dof[i]
            (x_dof[i] < 0) && (flow_direction = false)
        end

        if sym == :compressor
            ref(ss, sym, local_id)["flow"] = x_dof[i]
            (x_dof[i] < 0) && (compressor_flow_direction = false)
            ctrl_type, val = control(ss, :compressor, local_id)
            ref(ss, sym, local_id)["control_type"] = ctrl_type
            to_node = ref(ss, sym, local_id)["to_node"]
            fr_node = ref(ss, sym, local_id)["fr_node"]
            ref(ss, sym, local_id)["discharge_pressure"] =  x_dof[ref(ss, :node, to_node, :dof)]
            ref(ss, sym, local_id)["c_ratio"] = x_dof[ref(ss, :node, to_node, :dof)]/x_dof[ref(ss, :node, fr_node, :dof)]    
        end
    end

    return flow_direction, compressor_flow_direction
end


function populate_solution!(ss::Transienssimulator)
    sol = ss.sol
    units = params(ss, :units)

    function pressure_convertor(pu) 
        (units == 0) && (return pu * nominal_values(ss, :pressure)) 
        return pascal_to_psi(pu * nominal_values(ss, :pressure))
    end 

    function mass_flow_convertor(pu)
        kgps_to_mmscfd = get_kgps_to_mmscfd_conversion_factor(params(ss))
        (units == 0) && (return pu * nominal_values(ss, :mass_flow)) 
        return pu * nominal_values(ss, :mass_flow) * kgps_to_mmscfd
    end 
    
    for i in collect(keys(ref(ss, :node)))       
        sol["nodal_pressure"][i] = pressure_convertor(ref(ss, :node, i, "pressure"))
    end

    for i in collect(keys(ref(ss, :pipe)))
        sol["pipe_flow"][i] = mass_flow_convertor(ref(ss, :pipe, i, "mass_flow"))
    end
    
    for i in collect(keys(ref(ss, :compressor)))
        sol["compressor_flow"][i] = mass_flow_convertor(ref(ss, :compressor, i, "flow"))
    end

    return
end 