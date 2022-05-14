


"""function assembles the residuals"""
function assemble_residual!(ss::SteadySimulator, x_dof::AbstractArray, residual_edges::AbstractArray, residual_nonslack::AbstractArray, A_ns::SparseArrays.SparseMatrixCSC, val_ns::AbstractArray)
    num_ns =  ss.ref[:total_nonslack_vertices]
    num_edges = ss.ref[:total_edges]

    _eval_pipe_equations!(ss, x_dof, residual_edges)
    _eval_compressor_equations!(ss, x_dof, residual_edges)
    f_edges = x_dof[num_ns+1:end]
    residual_nonslack = A_ns * f_edges + val_ns
    return
end

"""function assembles the Jacobians"""
function assemble_mat(ss::SteadySimulator)
   
    num_ns =  ss.ref[:total_nonslack_vertices]
    num_edges = ss.ref[:total_edges]

    n =  2*length(ref(ss, :pipe)) + 2*length(ref(ss, :compressor)) # 2 entries for each edge
    i_vec, j_vec, k_vec, k_alpha_vec =  Vector{Int32}(), Vector{Int32}(), Vector{Float64}(), Vector{Float64}()
    sizehint!(i_vec, n)
    sizehint!(j_vec, n)
    sizehint!(k_vec, n)
    sizehint!(k_alpha_vec, n)


    # _eval_junction_equations_sparse_mat!(ss, x_dof, i_vec, j_vec, k_vec)
    _eval_pipe_equations_sparse_mat!(ss, i_vec, j_vec, k_vec, k_alpha_vec)
    _eval_compressor_equations_sparse_mat!(ss, i_vec, j_vec, k_vec, k_alpha_vec)

    # note we want ns x e edge incidence  which is transpose of vertex incidence hence jvec first
    A_ns = sparse(j_vec, i_vec, k_vec, num_ns, num_edges)
    A_alpha_ns = sparse(j_vec, i_vec, k_alpha_vec, num_ns, num_edges)

    return A_ns, A_alpha_ns
end

"""in place Jacobian computation for pipes"""
function _eval_pipe_equations_sparse_mat!(ss::SteadySimulator, i_vec::Array, j_vec::Array, k_vec::Array, k_alpha_vec::Array)
    @inbounds for (key, pipe) in ref(ss, :pipe)
        edge_dof = pipe["edge_dof"]  #shld be edge_id
        fr_node = pipe["fr_node"]  
        to_node = pipe["to_node"]


        if ref(ss, :node, fr_node, "vertex_dof") > 0
            push!(i_vec, edge_dof)
            push!(j_vec, ref(ss, :node, fr_node, "vertex_dof"))
            push!(k_vec, +1)
            push!(k_alpha_vec, +1)
        end

        if ref(ss, :node, to_node, "vertex_dof") > 0
            push!(i_vec, edge_dof)
            push!(j_vec, ref(ss, :node, to_node, "vertex_dof"))
            push!(k_vec, -1)
            push!(k_alpha_vec, -1)
        end

    end
end

"""in place Jacobian computation for compressors"""
function _eval_compressor_equations_sparse_mat!(ss::SteadySimulator, i_vec::Array, j_vec::Array, k_vec::Array, k_alpha_vec::Array)
    (!haskey(ref(ss), :compressor)) && (return)
    @inbounds for (comp_id, comp) in ref(ss, :compressor)
        edge_dof = comp["edge_dof"] 
        ctr, cmpr_val = control(ss, :compressor, comp_id) #assuming only cratio ctrl
        to_node = comp["to_node"]
        fr_node = comp["fr_node"]

        if ref(ss, :node, fr_node, "vertex_dof") > 0
            push!(i_vec, edge_dof)
            push!(j_vec, ref(ss, :node, fr_node, "vertex_dof"))
            push!(k_vec, +1)
            push!(k_alpha_vec, +1 * cmpr_val^2)
        end

        if ref(ss, :node, to_node, "vertex_dof") > 0
            push!(i_vec, edge_dof)
            push!(j_vec, ref(ss, :node, to_node, "vertex_dof"))
            push!(k_vec, -1)
            push!(k_alpha_vec, -1)
        end
        
    end
end



"""residual computation for pipes"""
function _eval_pipe_equations!(ss::SteadySimulator, x_dof::AbstractArray, residual_edges::AbstractArray)
    @inbounds for (_, pipe) in ref(ss, :pipe)
        dofid = pipe["dof"]
        edge_dof = pipe["edge_dof"]
        f = x_dof[dofid]
        fr_node = pipe["fr_node"]  
        to_node = pipe["to_node"]
        

        if ref(ss, :node, fr_node, "vertex_dof") > 0
            pi_fr = x_dof[ref(ss, :node, fr_node, "dof")]
        else
            ctrl_type, val = control(ss, :node, fr_node) # val is  pressure
            pi_fr = get_potential(ss, val)
        end

        if ref(ss, :node, to_node, "vertex_dof") > 0
            pi_to = x_dof[ref(ss, :node, to_node, "dof")]
        else
            ctrl_type, val = control(ss, :node, to_node) # val is  pressure
            pi_to = get_potential(ss, val)
        end

        c = nominal_values(ss, :mach_num)^2 / nominal_values(ss, :euler_num) 

        resistance = pipe["friction_factor"] * pipe["length"] * c / (2 * pipe["diameter"] * pipe["area"]^2)
        residual_edges[edge_dof] = pi_fr - pi_to - f * abs(f) * resistance
    end
end

"""residual computation for compressor"""
function _eval_compressor_equations!(ss::SteadySimulator, x_dof::AbstractArray, residual_edges::AbstractArray)
    (!haskey(ref(ss), :compressor)) && (return)
    @inbounds for (comp_id, comp) in ref(ss, :compressor)
        # dofid = comp["dof"] 
        edge_dof = comp["edge_dof"]

        ctr, cmpr_val = control(ss, :compressor, comp_id)
        to_node = comp["to_node"]
        fr_node = comp["fr_node"]

        if ref(ss, :node, fr_node, "vertex_dof") > 0
            pi_fr = x_dof[ref(ss, :node, fr_node, "dof")]
        else
            ctrl_type, val = control(ss, :node, fr_node) # val is  pressure
            pi_fr = get_potential(ss, val)
        end

        if ref(ss, :node, to_node, "vertex_dof") > 0
            pi_to = x_dof[ref(ss, :node, to_node, "dof")]
        else
            ctrl_type, val = control(ss, :node, to_node) # val is  pressure
            pi_to = get_potential(ss, val)
        end
        
            
        residual_edges[edge_dof] = (cmpr_val^2) * pi_fr -  pi_to
    end
end







