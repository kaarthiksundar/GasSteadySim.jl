function run_simulator!(ss::SteadySimulator; 
    method::Symbol=:newton,
    iteration_limit::Int64=2000, 
    kwargs...)
    

    ## Stage 1
    A_ns, A_alpha_ns, A_slack = assemble_mat(ss)
    println(A_ns, "\n", A_alpha_ns, "\n", A_slack)

    num_ns =  ss.ref[:total_nonslack_vertices]
    num_edges = ss.ref[:total_edges]
    # return A_ns, A_alpha_ns, num_ns, num_edges


    ## Stage 2
    x = Variable()
    xmax = 1
    y = Semidefinite(num_ns)
    z = Semidefinite(num_edges)
    # SDP constraints
    p = minimize(-x, (y + x) * A_ns - A_alpha_ns * (z + x)   == zeros(num_ns, num_edges), 
        x > 0,  x <= xmax, opnorm(y, 2) <= 10*xmax, opnorm(z, 2) <= 10*xmax )
    solve!(p, SCS.Optimizer; silent_solver = false)
    println(p.status)
    W1 = evaluate(y)
    W2= evaluate(z)


    
    # return p

    ## Stage 3
    ## next  block required to compute residual
    residual_nonslack = zeros(Float64, num_ns) 
    residual_edges = zeros(Float64, num_edges)
    val_ns = zeros(Float64, num_ns) 

    @inbounds for (node_id, junction) in ref(ss, :node)

        if junction["vertex_dof"] < 0
            continue
        end
        ctrl_type, val = control(ss, :node, node_id) # val is withdrawal 
        val_ns[junction["vertex_dof"]] = -val

    end

    #here do VI soln algo in loop
    x_k0 = ones(Float64, num_ns + num_edges)
    for i = 1: 5
        assemble_residual!(ss, x_k0, residual_edges, residual_nonslack, A_ns, val_ns, W1, W2)
        x_int[1:num_ns] = x_k0[1:num_ns] - residual_nonslack/1000 
        x_int[num_ns+1:num_ns+num_edges] = x_k0[num_ns:num_edges] - residual_edges/1000 
        #
        assemble_residual!(ss, x_int, residual_edges, residual_nonslack, A_ns, val_ns, W1, W2)
        x_k0[1:num_ns] = x_k0[1:num_ns] - residual_nonslack/1000 
        x_k0[num_ns+1:num_ns+num_edges] = x_k0[num_ns:num_edges] - residual_edges/1000 
    end


 


    

    # x_guess = _create_initial_guess_dof!(ss)
    # n = length(x_guess)
    

    # residual_fun! = (r_dof, x_dof) -> assemble_residual!(ss, x_dof, r_dof)
    # Jacobian_fun! = (J_dof, x_dof) -> assemble_mat!(ss, x_dof, J_dof)
    # J0 = spzeros(n, n)
    # assemble_mat!(ss, rand(n), J0)
    # df = OnceDifferentiable(residual_fun!, Jacobian_fun!, rand(n), rand(n), J0)

    # time = @elapsed soln = nlsolve(df, x_guess; method = method, iterations = iteration_limit, kwargs...)

    # convergence_state = converged(soln)

    # if convergence_state == false
    #     return SolverReturn(nl_solve_failure, 
    #         soln.iterations, 
    #         soln.residual_norm, 
    #         time, soln.zero, 
    #         Int[], Int[], Int[])
    # end

    # sol_return = update_solution_fields_in_ref!(ss, soln.zero)
    # populate_solution!(ss)

    # unphysical_solution = ~isempty(sol_return[:compressors_with_neg_flow]) || 
    # ~isempty(sol_return[:nodes_with_neg_potential])

    # if unphysical_solution
    #     is_unique = isempty(sol_return[:nodes_with_pressure_not_in_domain])
    #     if is_unique 
    #         return SolverReturn(unique_unphysical_solution, 
    #             soln.iterations, 
    #             soln.residual_norm, 
    #             time, soln.zero, 
    #             sol_return[:compressors_with_neg_flow], 
    #             sol_return[:nodes_with_neg_potential],
    #             sol_return[:nodes_with_pressure_not_in_domain])
    #     else 
    #         return SolverReturn(unphysical_solution, 
    #             soln.iterations, 
    #             soln.residual_norm, 
    #             time, soln.zero, 
    #             sol_return[:compressors_with_neg_flow], 
    #             sol_return[:nodes_with_neg_potential],
    #             sol_return[:nodes_with_pressure_not_in_domain])
    #     end 
    # end 

    
    # return SolverReturn(unique_physical_solution, 
    #     soln.iterations, 
    #     soln.residual_norm, 
    #     time, soln.zero, 
    #     sol_return[:compressors_with_neg_flow], 
    #     sol_return[:nodes_with_neg_potential],
    #     sol_return[:nodes_with_pressure_not_in_domain])

end

# function _create_initial_guess_dof!(ss::SteadySimulator)::Array
#     ndofs = length(ref(ss, :dof))
#     x_guess = 0.5 * ones(Float64, ndofs) 
#     dofs_updated = 0

#     components = [:node, :pipe, :compressor, 
#         :control_valve, :valve, 
#         :resistor, :loss_resistor, :short_pipe]

#     for component in components 
#         for (i, val) in get(ss.initial_guess, component, [])
#             x_guess[ref(ss, component, i, "dof")] = val 
#             dofs_updated += 1
#         end 
#     end 
#     return x_guess
# end