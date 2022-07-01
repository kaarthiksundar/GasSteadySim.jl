


"""function assembles the residuals"""
function assemble_residual!(ss::SteadySimulator, x_dof::AbstractArray, residual_edges::AbstractArray, residual_nonslack::AbstractArray, 
    A_ns::SparseArrays.SparseMatrixCSC, A_alpha_ns_transpose::SparseArrays.SparseMatrixCSC, A_slack::Array,
    val_ns::Array, beta_vec::Array, pi_slack_vec::Array, W1::AbstractArray, W2::AbstractArray)
    num_ns =  ss.new_ref[:total_nonslack_vertices]
    num_edges = ss.new_ref[:total_edges]

    pi_ns = x_dof[1:num_ns]
    f_edges = x_dof[num_ns+1:end]
    residual_nonslack .= W1 * (A_ns * f_edges + val_ns) #note the dot, else will not modify argument
    residual_edges .= W2 * (A_alpha_ns_transpose * pi_ns + transpose(A_slack) * pi_slack_vec - beta_vec .* f_edges .* abs.(f_edges))
    return
end

"""function assembles the Jacobians"""
function assemble_mat(ss::SteadySimulator)
   
    num_ns =  ss.new_ref[:total_nonslack_vertices]
    num_edges = ss.new_ref[:total_edges]
    num_slack = ss.new_ref[:total_slack_vertices]
    A_slack = zeros(num_slack, num_edges)
    beta_vec = zeros(num_edges)
    val_ns = zeros(num_ns)
    pi_slack_vec = zeros(num_slack)

    n =  2*length(new_ref(ss, :pipe))  # 2 entries for each edge
    i_vec, j_vec, k_vec, k_alpha_vec =  Vector{Int32}(), Vector{Int32}(), Vector{Float64}(), Vector{Float64}()
    sizehint!(i_vec, n)
    sizehint!(j_vec, n)
    sizehint!(k_vec, n)
    sizehint!(k_alpha_vec, n)


    _eval_pipe_equations_sparse_mat!(ss, i_vec, j_vec, k_vec, k_alpha_vec, A_slack, beta_vec, val_ns, pi_slack_vec)

    # note we want ns x e edge incidence  which is transpose of vertex incidence hence jvec first
    A_ns = sparse(j_vec, i_vec, k_vec, num_ns, num_edges)
    A_alpha_ns = sparse(j_vec, i_vec, k_alpha_vec, num_ns, num_edges)

    return A_ns, A_alpha_ns, A_slack, beta_vec, val_ns, pi_slack_vec
end

"""in place Jacobian computation for pipes"""
function _eval_pipe_equations_sparse_mat!(ss::SteadySimulator, i_vec::Array, j_vec::Array, k_vec::Array, 
    k_alpha_vec::Array, A_slack::Array, beta_vec::Array, val_ns::Array, pi_slack_vec::Array)
    @inbounds for (key, pipe) in new_ref(ss, :pipe)
        edge_dof = pipe["edge_dof"]  #shld be edge_id
        fr_node = pipe["fr_node"]  
        to_node = pipe["to_node"]
        fr_factor = pipe["fr_factor"]
        to_factor = pipe["to_factor"]

        c = nominal_values(ss, :mach_num)^2 / nominal_values(ss, :euler_num) 
        resistance = pipe["friction_factor"] * pipe["length"] * c / (2 * pipe["diameter"] * pipe["area"]^2)
        beta_vec[edge_dof] = resistance


        if new_ref(ss, :node, fr_node, "vertex_dof") > 0
            push!(i_vec, edge_dof)
            push!(j_vec, new_ref(ss, :node, fr_node, "vertex_dof"))
            push!(k_vec, +1)
            push!(k_alpha_vec, +fr_factor)
            val_ns[new_ref(ss, :node, fr_node, "vertex_dof")] = -1 * new_ref(ss, :node, fr_node, "withdrawal")
        else 
            A_slack[abs(new_ref(ss, :node, fr_node, "vertex_dof")), edge_dof] = +1
            pi_slack_vec[abs(new_ref(ss, :node, fr_node, "vertex_dof"))] = get_potential(ss, new_ref(ss, :node, fr_node, "pressure"))
        end

        if new_ref(ss, :node, to_node, "vertex_dof") > 0
            push!(i_vec, edge_dof)
            push!(j_vec, new_ref(ss, :node, to_node, "vertex_dof"))
            push!(k_vec, -1)
            push!(k_alpha_vec, -to_factor)
            val_ns[new_ref(ss, :node, to_node, "vertex_dof")] = -1 * new_ref(ss, :node, to_node, "withdrawal")
        else
            A_slack[abs(new_ref(ss, :node, to_node, "vertex_dof")), edge_dof] = -1
            pi_slack_vec[abs(new_ref(ss, :node, to_node, "vertex_dof"))] = get_potential(ss, new_ref(ss, :node, to_node, "pressure"))
        end

    end
    return
end









