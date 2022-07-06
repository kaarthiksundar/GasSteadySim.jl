


"""function assembles the residuals"""
function assemble_residual!(ss::SteadySimulator, x_dof::AbstractArray, residual_dof::AbstractArray; continuation_param::Float64=1.0)
    _eval_junction_equations!(ss, x_dof, residual_dof)
    _eval_pipe_equations!(ss, x_dof, residual_dof, continuation_param=continuation_param)
    _eval_compressor_equations!(ss, x_dof, residual_dof)
    _eval_control_valve_equations!(ss, x_dof, residual_dof)
    _eval_pass_through_equations!(ss, x_dof, residual_dof)
end

"""function assembles the Jacobians"""
function assemble_mat!(ss::SteadySimulator, x_dof::AbstractArray, J::AbstractArray; continuation_param::Float64=1.0)
    fill!(J, 0)
    _eval_junction_equations_mat!(ss, x_dof, J)
    _eval_pipe_equations_mat!(ss, x_dof, J, continuation_param=continuation_param)
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
            coeff = (ref(ss, :is_pressure_node, node_id)) ? val : get_potential(ss, val)
            residual_dof[eqn_no] = x_dof[eqn_no] - coeff

        end

        if  ctrl_type == flow_control
            r = (-val) # inflow is positive convention
            out_edge = ref(ss, :outgoing_dofs, node_id)
            in_edge = ref(ss, :incoming_dofs, node_id)
            r -= sum(x_dof[e] for e in out_edge; init=0.0) 
            r += sum(x_dof[e] for e in in_edge; init=0.0)
            residual_dof[eqn_no] = r
        end
    end
end

"""residual computation for pipes"""
function _eval_pipe_equations!(ss::SteadySimulator, x_dof::AbstractArray, residual_dof::AbstractArray; continuation_param::Float64=1.0)
    @inbounds for (_, pipe) in ref(ss, :pipe)
        eqn_no = pipe["dof"]
        f = x_dof[eqn_no]
        fr_node = pipe["fr_node"]  
        to_node = pipe["to_node"]
        fr_dof = ref(ss, :node, fr_node, "dof")
        to_dof = ref(ss, :node, to_node, "dof")
        is_fr_pressure_node = ref(ss, :is_pressure_node, fr_node)
        is_to_pressure_node = ref(ss, :is_pressure_node, to_node)
        pi_fr = (is_fr_pressure_node) ? get_potential(ss, x_dof[fr_dof]) : x_dof[fr_dof] 
        pi_to = (is_to_pressure_node) ? get_potential(ss, x_dof[to_dof]) : x_dof[to_dof] 
        c = nominal_values(ss, :mach_num)^2 / nominal_values(ss, :euler_num) 

        resistance = 1e-3 * pipe["friction_factor"] * pipe["length"] * c / (2 * pipe["diameter"] * pipe["area"]^2)
        residual_dof[eqn_no] = pi_fr - pi_to -  resistance * f * ( abs(f) ) ^ continuation_param
        # instead of f * abs(f), do f * ( abs(f) )^z where z can be varied from 0 to 1 
    end
end

"""residual computation for compressor"""
function _eval_compressor_equations!(ss::SteadySimulator, x_dof::AbstractArray, residual_dof::AbstractArray)
    (!haskey(ref(ss), :compressor)) && (return)
    @inbounds for (comp_id, comp) in ref(ss, :compressor)
        eqn_no = comp["dof"] 
        ctr, cmpr_val = control(ss, :compressor, comp_id)
        to_node = comp["to_node"]
        fr_node = comp["fr_node"]
        
        if ctr  == c_ratio_control
            is_pressure_eq = ref(ss, :is_pressure_node, fr_node) || ref(ss, :is_pressure_node, to_node)
            val = (is_pressure_eq) ? cmpr_val : cmpr_val^2
            residual_dof[eqn_no] = x_dof[ref(ss, :node, to_node, "dof")] - val * x_dof[ref(ss, :node, fr_node, "dof")]
        elseif ctr == flow_control
            residual_dof[eqn_no] = x_dof[eqn_no] - cmpr_val
        elseif ctr == discharge_pressure_control
            coeff = ref(ss, :is_pressure_node, to_node) ? cmpr_val : get_potential(ss, cmpr_val)
            residual_dof[eqn_no] = x_dof[ref(ss, :node, to_node, "dof")] - coeff 
        end
    end
end

"""residual computation for control_valves"""
function _eval_control_valve_equations!(ss::SteadySimulator, x_dof::AbstractArray, residual_dof::AbstractArray)
    (!haskey(ref(ss), :control_valve)) && (return)
    @inbounds for (cv_id, cv) in ref(ss, :control_valve)
        eqn_no = cv["dof"] 
        ctr, cv_val = control(ss, :control_valve, cv_id)
        to_node = cv["to_node"]
        fr_node = cv["fr_node"]
        
        if ctr  == c_ratio_control
            is_pressure_eq = ref(ss, :is_pressure_node, fr_node) || ref(ss, :is_pressure_node, to_node)
            val = (is_pressure_eq) ? cv_val : cv_val^2
            residual_dof[eqn_no] = x_dof[ref(ss, :node, to_node, "dof")] - val * x_dof[ref(ss, :node, fr_node, "dof")]
        elseif ctr == flow_control
            residual_dof[eqn_no] = x_dof[eqn_no] - cv_val
        elseif ctr == discharge_pressure_control
            coeff = ref(ss, :is_pressure_node, to_node) ? cv_val : get_potential(ss, cv_val)
            residual_dof[eqn_no] = x_dof[ref(ss, :node, to_node, "dof")] - coeff
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
            fr_dof = ref(ss, :node, fr_node, "dof")
            to_dof = ref(ss, :node, to_node, "dof")
            is_fr_pressure_node = ref(ss, :is_pressure_node, fr_node)
            is_to_pressure_node = ref(ss, :is_pressure_node, to_node)
            if (is_fr_pressure_node && is_to_pressure_node)
                residual_dof[eqn_no] = x_dof[to_dof] - x_dof[fr_dof]
            else
                pi_fr = (is_fr_pressure_node) ? get_potential(ss, x_dof[fr_dof]) : x_dof[fr_dof] 
                pi_to = (is_to_pressure_node) ? get_potential(ss, x_dof[to_dof]) : x_dof[to_dof] 
                residual_dof[eqn_no] = pi_to - pi_fr
            end 
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
        J::AbstractArray; continuation_param::Float64=1.0)
    @inbounds for (key, pipe) in ref(ss, :pipe)
        eqn_no = pipe["dof"] 
        f = x_dof[eqn_no]
        fr_node = pipe["fr_node"]  
        to_node = pipe["to_node"]

        eqn_fr = ref(ss, :node, fr_node, "dof")
        eqn_to = ref(ss, :node, to_node, "dof")
        is_fr_pressure_node = ref(ss, :is_pressure_node, fr_node)
        is_to_pressure_node = ref(ss, :is_pressure_node, to_node)
        
        c = nominal_values(ss, :mach_num)^2 / nominal_values(ss, :euler_num) 
        resistance = 1e-3 * pipe["friction_factor"] * pipe["length"] * c / (2 * pipe["diameter"] * pipe["area"]^2)

        pi_dash_fr = (is_fr_pressure_node) ? get_potential_derivative(ss, x_dof[eqn_fr]) : 1.0 
        pi_dash_to = (is_to_pressure_node) ? get_potential_derivative(ss, x_dof[eqn_to]) : 1.0 

        J[eqn_no, eqn_fr] = pi_dash_fr
        J[eqn_no, eqn_to] = -pi_dash_to
        J[eqn_no, eqn_no] = -(1 + continuation_param) * resistance * ( abs(f) ) ^ continuation_param
        # use -(1+z) * ( abs(f) )^z so that for z = 1 we get -2 * abs(f)
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
        is_pressure_eq = ref(ss, :is_pressure_node, fr_node) || ref(ss, :is_pressure_node, to_node)
        
        if ctr  == c_ratio_control
            J[eqn_no, eqn_to] = 1
            J[eqn_no, eqn_fr] = (is_pressure_eq) ? (-cmpr_val) : (-cmpr_val^2)
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
        is_pressure_eq = ref(ss, :is_pressure_node, fr_node) || ref(ss, :is_pressure_node, to_node)
        
        if ctr  == c_ratio_control
            J[eqn_no, eqn_to] = 1
            J[eqn_no, eqn_fr] = is_pressure_eq ? (-cv_val) : (-cv_val^2)
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
            is_fr_pressure_node = ref(ss, :is_pressure_node, fr_node)
            is_to_pressure_node = ref(ss, :is_pressure_node, to_node)

            if (is_fr_pressure_node && is_to_pressure_node)
                J[eqn_no, eqn_to] = 1.0
                J[eqn_no, eqn_fr] = -1.0
            else
                pi_dash_fr = (is_fr_pressure_node) ? get_potential_derivative(ss, x_dof[eqn_fr]) : 1.0 
                pi_dash_to = (is_to_pressure_node) ? get_potential_derivative(ss, x_dof[eqn_to]) : 1.0
                J[eqn_no, eqn_to] = pi_dash_to
                J[eqn_no, eqn_fr] = -pi_dash_fr
            end 
        end 
    end 
end 