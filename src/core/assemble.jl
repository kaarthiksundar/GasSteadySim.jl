


"""function assembles the residuals"""
function assemble_residual!(ss::SteadySimulator, x_dof::AbstractArray, residual_dof::AbstractArray)
    _eval_junction_equations!(ss, x_dof, residual_dof)
    _eval_pipe_equations!(ss, x_dof, residual_dof)
    _eval_compressor_equations!(ss, x_dof, residual_dof)
    _eval_control_valve_equations!(ss, x_dof, residual_dof)
    _eval_pass_through_equations!(ss, x_dof, residual_dof)
end

"""function assembles the Jacobians"""
function assemble_mat!(ss::SteadySimulator, x_dof::AbstractArray, J::AbstractArray)
    fill!(J, 0)
    _eval_junction_equations_mat!(ss, x_dof, J)
    _eval_pipe_equations_mat!(ss, x_dof, J)
    _eval_compressor_equations_mat!(ss, x_dof, J)
    _eval_control_valve_equations_mat!(ss, x_dof, J)
    _eval_pass_through_equations_mat!(ss, x_dof, J)
end

"""residual computation for junctions"""
function _eval_junction_equations!(ss::SteadySimulator, x_dof::AbstractArray, residual_dof::AbstractArray)
    @inbounds for (node_id, junction) in ref(ss, :node)
        eqn_no = junction["dof"]
        ctrl_type, val = control(ss, :node, node_id) # val is withdrawal or pressure

        if ctrl_type == pressure_control
            residual_dof[eqn_no] = x_dof[eqn_no] - val
        end

        if  ctrl_type == flow_control
            r = (-val) # inflow is positive convention
            out_edge = ref(ss, :outgoing_dofs, node_id)
            in_edge = ref(ss, :incoming_dofs, node_id)
            r -= sum(x_dof[e] for e in out_edge) 
            r += sum(x_dof[e] for e in in_edge)
            residual_dof[eqn_no] = r
        end
    end
end

"""residual computation for pipes"""
function _eval_pipe_equations!(ss::SteadySimulator, x_dof::AbstractArray, residual_dof::AbstractArray)
    @inbounds for (_, pipe) in ref(ss, :pipe)
        eqn_no = pipe["dof"]
        f = x_dof[eqn_no]
        fr_node = pipe["fr_node"]  
        to_node = pipe["to_node"]
        c = nominal_values(ss, :mach_num)^2 / nominal_values(ss, :euler_num) 

        b1, b2 = get_eos_coeffs(ss)
        pressure_sqr_diff = x_dof[ref(ss, :node, fr_node, "dof")]^2 - x_dof[ref(ss, :node, to_node, "dof")]^2
        pressure_cube_diff = x_dof[ref(ss, :node, fr_node, "dof")]^3 - x_dof[ref(ss, :node, to_node, "dof")]^3 
        resistance = pipe["friction_factor"] * pipe["length"] * c / (2 * pipe["diameter"] * pipe["area"]^2)
        residual_dof[eqn_no] = (b1/2) * pressure_sqr_diff + (b2/3) * pressure_cube_diff - f * abs(f) * resistance
    end
end

"""residual computation for compressor"""
function _eval_compressor_equations!(ss::SteadySimulator, x_dof::AbstractArray, residual_dof::AbstractArray)
    (!haskey(ref(ss), :compressor)) && (return)
    @inbounds for (comp_id, comp) in ref(ss, :compressor)
        eqn_no = comp["dof"] 
        ctr, cmpr_val = control(ss, :compressor, comp_id)
        
        if ctr  == c_ratio_control
            to_node = comp["to_node"]
            fr_node = comp["fr_node"]
            residual_dof[eqn_no] = x_dof[ref(ss, :node, to_node, "dof")] - cmpr_val * x_dof[ref(ss, :node, fr_node,"dof")]
        elseif ctr == flow_control
            residual_dof[eqn_no] = x_dof[eqn_no] - cmpr_val
        elseif ctr == discharge_pressure_control
            to_node = comp["to_node"]
            residual_dof[eqn_no] = x_dof[ref(ss, :node, to_node, "dof")] - cmpr_val
        end
    end
end

"""residual computation for control_valves"""
function _eval_control_valve_equations!(ss::SteadySimulator, x_dof::AbstractArray, residual_dof::AbstractArray)
    (!haskey(ref(ss), :control_valve)) && (return)
    @inbounds for (cv_id, cv) in ref(ss, :control_valve)
        eqn_no = cv["dof"] 
        ctr, cv_val = control(ss, :control_valve, cv_id)
        
        if ctr  == c_ratio_control
            to_node = cv["to_node"]
            fr_node = cv["fr_node"]
            residual_dof[eqn_no] = x_dof[ref(ss, :node, to_node, "dof")] - 
                cv_val * x_dof[ref(ss, :node, fr_node,"dof")]
        elseif ctr == flow_control
            residual_dof[eqn_no] = x_dof[eqn_no] - cv_val
        elseif ctr == discharge_pressure_control
            to_node = comp["to_node"]
            residual_dof[eqn_no] = x_dof[ref(ss, :node, to_node, "dof")] - cv_val
        end
    end
end

"""residual computation for pass through components"""
function _eval_pass_through_equations!(ss::SteadySimulator, x_dof::AbstractArray, residual_dof::AbstractArray)
    components = [:valve, :resistor, :loss_resistor, :short_pipe]
    @inbounds for component in components 
        (!haskey(ref(ss), component)) && (continue)
        for (_, comp) in ref(ss, component)
            eqn_no = comp["dof"]
            fr_node = comp["fr_node"]
            to_node = comp["to_node"]
            residual_dof[eqn_no] = x_dof[ref(ss, :node, to_node, "dof")] - x_dof[ref(ss, :node, fr_node, "dof")]
        end 
    end 
end 

"""in place Jacobian computation for junctions"""
function _eval_junction_equations_mat!(ss::SteadySimulator, x_dof::AbstractArray, 
        J::AbstractArray)
    @inbounds for (node_id, junction) in ref(ss, :node)
        eqn_no = junction["dof"]
        ctrl_type, _ = control(ss, :node, node_id) # val is withdrawal or pressure

        if ctrl_type == pressure_control
            J[eqn_no, eqn_no] = 1
            continue
        end

        if  ctrl_type == flow_control
            out_edge = ref(ss, :outgoing_dofs, node_id)
            in_edge = ref(ss, :incoming_dofs, node_id)
            for e in out_edge
                J[eqn_no, e] = -1
            end
            for e in in_edge
                J[eqn_no, e] = +1
            end
        end
    end
end


"""in place Jacobian computation for pipes"""
function _eval_pipe_equations_mat!(ss::SteadySimulator, x_dof::AbstractArray, 
        J::AbstractArray)
    @inbounds for (key, pipe) in ref(ss, :pipe)
        eqn_no = pipe["dof"] 
        f = x_dof[eqn_no]
        fr_node = pipe["fr_node"]  
        to_node = pipe["to_node"]
        eqn_to = ref(ss, :node, to_node, "dof")
        eqn_fr = ref(ss, :node, fr_node, "dof")
        p_fr = x_dof[eqn_fr]
        p_to = x_dof[eqn_to]
        c = nominal_values(ss, :mach_num)^2 / nominal_values(ss, :euler_num) 


        b1, b2 = get_eos_coeffs(ss)
        resistance = pipe["friction_factor"] * pipe["length"] * c / (2 * pipe["diameter"] * pipe["area"]^2)
        var1 = (b1/2) * 2 * p_fr + (b2/3) * 3 * p_to^2
        var2 = - (b1/2)* 2 * p_to - (b2/3) * 3 * p_to^2
        var3 =  -2.0 * f * sign(f) * resistance
        
        J[eqn_no, eqn_fr] = var1
        J[eqn_no, eqn_to] = var2
        J[eqn_no, eqn_no] = var3 
    end
end

"""in place Jacobian computation for compressors"""
function _eval_compressor_equations_mat!(ss::SteadySimulator, x_dof::AbstractArray, 
        J::AbstractArray)
    (!haskey(ref(ss), :compressor)) && (return)
    @inbounds for (comp_id, comp) in ref(ss, :compressor)
        eqn_no = comp["dof"] 
        ctr, cmpr_val = control(ss, :compressor, comp_id)
        to_node = comp["to_node"]
        fr_node = comp["fr_node"]
        eqn_to = ref(ss, :node, to_node, "dof")
        eqn_fr = ref(ss, :node, fr_node, "dof")
        
        if ctr  == c_ratio_control
            J[eqn_no, eqn_to] = 1
            J[eqn_no, eqn_fr] = (-cmpr_val)
        elseif ctr == flow_control
            J[eqn_no, eqn_no] = 1
        elseif ctr == discharge_pressure_control
            J[eqn_no, eqn_to] = 1
        end
    end
end

"""in place Jacobian computation for control_valves"""
function _eval_control_valve_equations_mat!(ss::SteadySimulator, x_dof::AbstractArray, 
        J::AbstractArray)
    (!haskey(ref(ss), :control_valve)) && (return)
    @inbounds for (cv_id, cv) in ref(ss, :control_valve)
        eqn_no = cv["dof"] 
        ctr, cv_val = control(ss, :control_valve, cv_id)
        to_node = cv["to_node"]
        fr_node = cv["fr_node"]
        eqn_to = ref(ss, :node, to_node, "dof")
        eqn_fr = ref(ss, :node, fr_node, "dof")
        
        if ctr  == c_ratio_control
            J[eqn_no, eqn_to] = 1
            J[eqn_no, eqn_fr] = (-cv_val)
        elseif ctr == flow_control
            J[eqn_no, eqn_no] = 1
        elseif ctr == discharge_pressure_control
            J[eqn_no, eqn_to] = 1
        end
    end
end

"""in place Jacobian computation for pass through components"""
function _eval_pass_through_equations_mat!(ss::SteadySimulator, x_dof::AbstractArray, 
        J::AbstractArray)
    components = [:valve, :resistor, :loss_resistor, :short_pipe]
    @inbounds for component in components 
        (!haskey(ref(ss), component)) && (continue)
        for (_, comp) in ref(ss, component)
            eqn_no = comp["dof"]
            to_node = comp["to_node"]
            fr_node = comp["fr_node"]
            eqn_to = ref(ss, :node, to_node, "dof")
            eqn_fr = ref(ss, :node, fr_node, "dof")

            J[eqn_no, eqn_to] = 1 
            J[eqn_no, eqn_fr] = 1.0 
        end 
    end 
end 