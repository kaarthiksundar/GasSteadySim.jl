


"""function assembles the residuals for in place computation"""
function assemble_residual_in_place(ss::SteadySimulator, x_dof::Array, residual_dof::Array)
    _eval_junction_equations!(ss, x_dof, residual_dof)
    _eval_pipe_equations!(ss, x_dof, residual_dof)
    _eval_compressor_equations!(ss, x_dof, residual_dof)
end

"""function assembles the Jacobians for in place computation"""
function assemble_mat_in_place(ss::SteadySimulator, x_dof::Array, J::Array)
    fill!(J, 0)
    _eval_junction_equations_mat!(ss, x_dof, J)
    _eval_pipe_equations_mat!(ss, x_dof, J)
    _eval_compressor_equations_mat!(ss, x_dof, J)
end

"""function assembles the residuals for sparse computation"""
function assemble_residual(ss::SteadySimulator, x_dof::Array)

    num_dofs = length(ref(ss, :dof))
    residual_dof = zeros(Float64, num_dofs)
    _eval_junction_equations!(ss, x_dof, residual_dof)
    _eval_pipe_equations!(ss, x_dof, residual_dof)
    _eval_compressor_equations!(ss, x_dof, residual_dof)
    return residual_dof
end

"""function assembles the Jacobians for sparse computation"""
function assemble_mat(ss::SteadySimulator, x_dof::Array)
    num_dofs = length(ref(ss, :dof))
    n =  length(ref(ss, :node)) + 3*length(ref(ss, :pipe)) + 3*length(ref(ss, :compressor))
    i_vec, j_vec, k_vec =  Vector{Int32}(), Vector{Int32}(), Vector{Float64}()
    sizehint!(i_vec, n)
    sizehint!(j_vec, n)
    sizehint!(k_vec, n)

    _eval_junction_equations_sparse_mat!(ss, x_dof, i_vec, j_vec, k_vec)
    _eval_pipe_equations_sparse_mat!(ss, x_dof, i_vec, j_vec, k_vec)
    _eval_compressor_equations_sparse_mat!(ss, x_dof, i_vec, j_vec, k_vec)
    J = sparse(i_vec, j_vec, k_vec, num_dofs, num_dofs)
    return J
end

"""residual computation for junctions"""
function _eval_junction_equations!(ss::SteadySimulator, x_dof::Array, residual_dof::Array)
    for (node_id, junction) in ref(ss, :node)
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
function _eval_pipe_equations!(ss::SteadySimulator, x_dof::Array, residual_dof::Array)
    for (_, pipe) in ref(ss, :pipe)
        eqn_no = pipe[:dof] 
        f = x_dof[eqn_no]
        fr_node = pipe["fr_node"]  
        to_node = pipe["to_node"]
        c = ss.nominal_values[:mass_flow]^2  / (ss.nominal_values[:pressure] * ss.nominal_values[:density])

        b1, b2 = get_eos_coeffs(ss)
        pressure_sqr_diff = x_dof[ref(ss, :node, fr_node, :dof)]^2 - x_dof[ref(ss, :node, to_node, :dof)]^2
        pressure_cube_diff = x_dof[ref(ss, :node, fr_node,:dof)]^3 - x_dof[ref(ss, :node, to_node,:dof)]^3
        resistance = pipe["friction_factor"] * pipe["length"] * c / (2 * pipe["diameter"] * pipe["area"]^2)
        residual_dof[eqn_no] = (b1/2) * pressure_sqr_diff + (b2/3) * pressure_cube_diff - f * abs(f) * resistance
    end
end

"""residual computation for compressor"""
function _eval_compressor_equations!(ss::SteadySimulator, x_dof::Array, residual_dof::Array)
    for (comp_id, comp) in ref(ss, :compressor)
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
function _eval_junction_equations_mat!(ss::SteadySimulator, x_dof::Array, J::Array)
    for (node_id, junction) in ref(ss, :node)
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
function _eval_pipe_equations_mat!(ss::SteadySimulator, x_dof::Array, J::Array)
    for (key, pipe) in ref(ss, :pipe)
        eqn_no = pipe[:dof] 
        f = x_dof[eqn_no]
        fr_node = pipe["fr_node"]  
        to_node = pipe["to_node"]
        eqn_to = ref(ss, :node, to_node, :dof)
        eqn_fr = ref(ss, :node, fr_node, :dof)
        p_fr = x_dof[eqn_fr]
        p_to = x_dof[eqn_to]
        c = ss.nominal_values[:mass_flow]^2  / (ss.nominal_values[:pressure] * ss.nominal_values[:density])

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
function _eval_compressor_equations_mat!(ss::SteadySimulator, x_dof::Array, J::Array)
    for (comp_id, comp) in ref(ss, :compressor)
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

"""sparse Jacobian computation for junctions"""
function _eval_junction_equations_sparse_mat!(ss::SteadySimulator, x_dof::Array, i_vec::Array, j_vec::Array, k_vec::Array)
    for (node_id, junction) in ref(ss, :node)
        eqn_no = junction[:dof]
        ctrl_type, _ = control(ss, :node, node_id) #val is withdrawal or pressure

        if ctrl_type == pressure_control
            push!(i_vec, eqn_no)
            push!(j_vec, eqn_no)
            push!(k_vec, 1)
        end

        if  ctrl_type == flow_control
            out_edge = ref(ss, :outgoing_dofs, node_id)
            in_edge = ref(ss, :incoming_dofs, node_id)
            for e in out_edge
                push!(i_vec, eqn_no)
                push!(j_vec, e)
                push!(k_vec, -1)
            end
            for e in in_edge
                push!(i_vec, eqn_no)
                push!(j_vec, e)
                push!(k_vec, +1)
            end
        end
    end
end

"""sparse Jacobian computation for pipes"""   
function _eval_pipe_equations_sparse_mat!(ss::SteadySimulator, x_dof::Array, i_vec::Array, j_vec::Array, k_vec::Array)
    for (key, pipe) in ref(ss, :pipe)
        eqn_no = pipe[:dof] 
        fr_node = pipe["fr_node"]  
        to_node = pipe["to_node"]
        eqn_to = ref(ss, :node, to_node, :dof)
        eqn_fr = ref(ss, :node, fr_node, :dof)
        f = x_dof[eqn_no]
        p_fr = x_dof[eqn_fr]
        p_to = x_dof[eqn_to]
        c = ss.nominal_values[:mass_flow]^2  / (ss.nominal_values[:pressure] * ss.nominal_values[:density])

        b1, b2 = get_eos_coeffs(ss)
        resistance = pipe["friction_factor"] * pipe["length"] * c / (2 * pipe["diameter"] * pipe["area"]^2)
        var1 = b1 * p_fr + b2 * p_fr^2
        var2 = -b1 * p_to - b2 * p_to^2
        var3 =  - 2 * f * sign(f) * resistance

        push!(i_vec, eqn_no)
        push!(j_vec, eqn_fr)
        push!(k_vec, var1)

        push!(i_vec, eqn_no)
        push!(j_vec, eqn_to)
        push!(k_vec, var2)
        
        push!(i_vec, eqn_no)
        push!(j_vec, eqn_no)
        push!(k_vec, var3)
    end
end

"""sparse Jacobian computation for compressors""" 
function _eval_compressor_equations_sparse_mat!(ss::SteadySimulator, x_dof::Array, i_vec::Array, j_vec::Array, k_vec::Array)
    
    for (comp_id, comp) in ref(ss, :compressor)
        eqn_no = comp[:dof] 
        to_node = comp["to_node"]
        fr_node = comp["fr_node"]
        eqn_to = ref(ss, :node, to_node, :dof)
        eqn_fr = ref(ss, :node, fr_node, :dof)
        
        ctr, cmpr_val = control(ss, :compressor, comp_id)

        if ctr  == c_ratio_control
            push!(i_vec, eqn_no)
            push!(j_vec, eqn_to)
            push!(k_vec, 1)
            push!(i_vec, eqn_no)
            push!(j_vec, eqn_fr)
            push!(k_vec, -cmpr_val)
        elseif ctr == flow_control
            push!(i_vec, eqn_no)
            push!(j_vec, eqn_no)
            push!(k_vec, 1)
        elseif ctr == discharge_pressure_control
            push!(i_vec, eqn_no)
            push!(j_vec, eqn_to)
            push!(k_vec, 1)
        end
    end
end




