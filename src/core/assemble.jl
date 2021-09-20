


"""function assembles the residuals"""
function assemble_residual!(ss::SteadySimulator, x_dof::AbstractArray, residual_dof::AbstractArray)
    _eval_junction_equations!(ss, x_dof, residual_dof)
    _eval_pipe_equations!(ss, x_dof, residual_dof)
    _eval_compressor_equations!(ss, x_dof, residual_dof)
end

"""function assembles the Jacobians"""
function assemble_mat!(ss::SteadySimulator, x_dof::AbstractArray, J::AbstractArray)
    fill!(J, 0)
    _eval_junction_equations_mat!(ss, x_dof, J)
    _eval_pipe_equations_mat!(ss, x_dof, J)
    _eval_compressor_equations_mat!(ss, x_dof, J)
end

"""residual computation for junctions"""
function _eval_junction_equations!(ss::SteadySimulator, x_dof::AbstractArray, residual_dof::AbstractArray)
    @inbounds for (node_id, junction) in ref(ss, :node)
        eqn_no = junction[:dof]
        ctrl_type, val = control(ss, :node, node_id) # val is withdrawal or pressure

        if ctrl_type == pressure_control
            residual_dof[eqn_no] = x_dof[eqn_no] - val
        end

        if  ctrl_type == flow_control
            r = (-val) # inflow is positive convention
            out_edge = ref(ss, :outgoing_dofs, node_id)
            in_edge = ref(ss, :incoming_dofs, node_id)
            for e in out_edge
                r -= x_dof[e]
            end
            for e in in_edge
                r += x_dof[e]
            end
            residual_dof[eqn_no] = r
        end
    end
end

"""residual computation for pipes"""
function _eval_pipe_equations!(ss::SteadySimulator, x_dof::AbstractArray, residual_dof::AbstractArray)
    @inbounds for (_, pipe) in ref(ss, :pipe)
        eqn_no = pipe[:dof] 
        f = x_dof[eqn_no]
        fr_node = pipe["fr_node"]  
        to_node = pipe["to_node"]
        c = nominal_values(ss, :mach_num)^2 / nominal_values(ss, :euler_num) 

        b1, b2 = get_eos_coeffs(ss)
        pressure_sqr_diff = x_dof[ref(ss, :node, fr_node, :dof)]^2 - x_dof[ref(ss, :node, to_node, :dof)]^2
        pressure_cube_diff = x_dof[ref(ss, :node, fr_node,:dof)]^3 - x_dof[ref(ss, :node, to_node,:dof)]^3
        resistance = pipe["friction_factor"] * pipe["length"] * c / (2 * pipe["diameter"] * pipe["area"]^2)
        residual_dof[eqn_no] = (b1/2) * pressure_sqr_diff + (b2/3) * pressure_cube_diff - f * abs(f) * resistance
    end
end

"""residual computation for compressor"""
function _eval_compressor_equations!(ss::SteadySimulator, x_dof::AbstractArray, residual_dof::AbstractArray)
    @inbounds for (comp_id, comp) in ref(ss, :compressor)
        eqn_no = comp[:dof] 
        ctr, cmpr_val = control(ss, :compressor, comp_id)
        
        if ctr  == c_ratio_control
            to_node = comp["to_node"]
            fr_node = comp["fr_node"]
            residual_dof[eqn_no] = x_dof[ref(ss, :node, to_node, :dof)] - cmpr_val * x_dof[ref(ss, :node, fr_node,:dof)]
        elseif ctr == flow_control
            residual_dof[eqn_no] = x_dof[eqn_no] - cmpr_val
        elseif ctr == discharge_pressure_control
            to_node = comp["to_node"]
            residual_dof[eqn_no] = x_dof[ref(ss, :node, to_node,:dof)] - cmpr_val
        end
    end
end


"""in place Jacobian computation for junctions"""
function _eval_junction_equations_mat!(ss::SteadySimulator, x_dof::AbstractArray, 
        J::AbstractArray)
    @inbounds for (node_id, junction) in ref(ss, :node)
        eqn_no = junction[:dof]
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
        eqn_no = pipe[:dof] 
        f = x_dof[eqn_no]
        fr_node = pipe["fr_node"]  
        to_node = pipe["to_node"]
        eqn_to = ref(ss, :node, to_node, :dof)
        eqn_fr = ref(ss, :node, fr_node, :dof)
        p_fr = x_dof[eqn_fr]
        p_to = x_dof[eqn_to]
        c = nominal_values(ss, :mach_num)^2 / nominal_values(ss, :euler_num) 


        b1, b2 = get_eos_coeffs(ss)
        resistance = pipe["friction_factor"] * pipe["length"] * c / (2 * pipe["diameter"] * pipe["area"]^2)
        var1 = b1 * p_fr + b2 * p_fr^2
        var2 = -b1 * p_to - b2 * p_to^2
        var3 =  -2.0 * f * sign(f) * resistance
        
        J[eqn_no, eqn_fr] = var1
        J[eqn_no, eqn_to] = var2
        J[eqn_no, eqn_no] = var3 
    end
end

"""in place Jacobian computation for compressors"""
function _eval_compressor_equations_mat!(ss::SteadySimulator, x_dof::AbstractArray, 
        J::AbstractArray)
    @inbounds for (comp_id, comp) in ref(ss, :compressor)
        eqn_no = comp[:dof] 
        ctr, cmpr_val = control(ss, :compressor, comp_id)
        to_node = comp["to_node"]
        fr_node = comp["fr_node"]
        eqn_to = ref(ss, :node, to_node, :dof)
        eqn_fr = ref(ss, :node, fr_node, :dof)
        
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