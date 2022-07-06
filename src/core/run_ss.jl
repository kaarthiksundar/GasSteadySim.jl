function run_simulator!(ss::SteadySimulator; 
    continuation_param::Float64=1.0,
    method::Symbol=:newton,
    iteration_limit::Int64=2000, 
    kwargs...)::SolverReturn
    
    x_guess = _create_initial_guess_dof!(ss)
    n = length(x_guess)

    residual_fun! = (r_dof, x_dof) -> assemble_residual!(ss, x_dof, r_dof, continuation_param=continuation_param)
    Jacobian_fun! = (J_dof, x_dof) -> assemble_mat!(ss, x_dof, J_dof, continuation_param=continuation_param)
    J0 = spzeros(n, n)
    assemble_mat!(ss, rand(n), J0, continuation_param=continuation_param)
    df = OnceDifferentiable(residual_fun!, Jacobian_fun!, rand(n), rand(n), J0)

    time = @elapsed soln = nlsolve(df, x_guess; method = method, iterations = iteration_limit, kwargs...)

    convergence_state = converged(soln)

    if convergence_state == false
        return SolverReturn(nl_solve_failure, 
            soln.iterations, 
            soln.residual_norm, 
            time, soln.zero, 
            Int[], Int[], Int[])
    end

    writedlm("dummy_initial.txt", soln.zero)

    sol_return = update_solution_fields_in_ref!(ss, soln.zero)
    populate_solution!(ss)

    unphysical_solution = ~isempty(sol_return[:compressors_with_neg_flow]) || 
    ~isempty(sol_return[:nodes_with_neg_potential])

    if unphysical_solution
        is_unique = isempty(sol_return[:nodes_with_pressure_not_in_domain])
        if is_unique 
            return SolverReturn(unique_unphysical_solution, 
                soln.iterations, 
                soln.residual_norm, 
                time, soln.zero, 
                sol_return[:compressors_with_neg_flow], 
                sol_return[:nodes_with_neg_potential],
                sol_return[:nodes_with_pressure_not_in_domain])
        else 
            return SolverReturn(unphysical_solution, 
                soln.iterations, 
                soln.residual_norm, 
                time, soln.zero, 
                sol_return[:compressors_with_neg_flow], 
                sol_return[:nodes_with_neg_potential],
                sol_return[:nodes_with_pressure_not_in_domain])
        end 
    end 

    
    return SolverReturn(unique_physical_solution, 
        soln.iterations, 
        soln.residual_norm, 
        time, soln.zero, 
        sol_return[:compressors_with_neg_flow], 
        sol_return[:nodes_with_neg_potential],
        sol_return[:nodes_with_pressure_not_in_domain])
end

function _create_initial_guess_dof!(ss::SteadySimulator)::Array
    ndofs = length(ref(ss, :dof))
     
    x_guess = readdlm("dummy_initial.txt")
    if length(x_guess) != ndofs
        x_guess = ones(Float64, ndofs)
    end

    # x_guess = rand(Normal(0, 1), ndofs)

    # num_nodes = length(ss.ref[:node])
    # x_guess[1:num_nodes] .= 1.0

    dofs_updated = 0

    components = [:node, :pipe, :compressor, 
        :control_valve, :valve, 
        :resistor, :loss_resistor, :short_pipe]

    for component in components 
        for (i, val) in get(ss.initial_guess, component, [])
            x_guess[ref(ss, component, i, "dof")] = val 
            dofs_updated += 1
        end 
    end 
    return x_guess
end