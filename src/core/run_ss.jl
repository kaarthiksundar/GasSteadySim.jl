function is_feasible!(ss::SteadySimulator, optimizer)::Bool

    set_optimizer(ss.feasibility_model, optimizer)
    optimize!(ss.feasibility_model)
    is_infeasible = (termination_status(ss.feasibility_model) == MOI.INFEASIBLE)
    return !is_infeasible
end

function run_simulator!(ss::SteadySimulator; 
    method::Symbol=:trust_region,
    iteration_limit::Int64=2000)::SolverReturn
    
    x_guess = _create_initial_guess_dof!(ss)
    n = length(x_guess)
    

    residual_fun! = (r_dof, x_dof) -> assemble_residual!(ss, x_dof, r_dof)
    Jacobian_fun! = (J_dof, x_dof) -> assemble_mat!(ss, x_dof, J_dof)
    J0 = spzeros(n, n)
    assemble_mat!(ss, rand(n), J0)
    df = OnceDifferentiable(residual_fun!, Jacobian_fun!, rand(n), rand(n), J0)

    t_first = @elapsed soln = nlsolve(df, x_guess; method = method, iterations = iteration_limit)

    t_second = 0.0
    all_pressures_non_neg = check_for_negative_pressures(ss, soln.zero)
    convergence_state = converged(soln)

    if convergence_state == false
        return SolverReturn(initial_nl_solve_failure, 
            soln.iterations, 
            soln.residual_norm, 
            t_first + t_second, 
            soln.zero, Int[])
    end

    if all_pressures_non_neg == true
        @info "initial nonlinear solve complete and all pressures non-negative"
    end

    iters_initial = soln.iterations

    if all_pressures_non_neg == false
        @info "correcting pressures..."
        reinitialize_for_positive_pressure!(ss, soln.zero)
        t_second = @elapsed soln = nlsolve(df, soln.zero; method = method, iterations = iteration_limit)
        
        all_pressures_non_neg = check_for_negative_pressures(ss, soln.zero)
        convergence_state = converged(soln)

        if convergence_state == false
            @warn "nonlinear solve for pressure correction failed!"
            update_solution_fields_in_ref!(ss, soln.zero)
            populate_solution!(ss)
            return SolverReturn(pressure_correction_nl_solve_failure, 
                iters_initial + soln.iterations, 
                soln.residual_norm, 
                t_first + t_second, 
                soln.zero, Int[])
        end

        if all_pressures_non_neg == false
            @warn "pressure correction failed"
            update_solution_fields_in_ref!(ss, soln.zero)
            populate_solution!(ss)
            return SolverReturn(pressure_correction_failure, 
                iters_initial + soln.iterations, 
                soln.residual_norm, 
                t_first + t_second, 
                soln.zero, Int[])
        else 
            @info "pressure correction successful"
        end
    end


    flow_direction, negative_flow_in_compressors = update_solution_fields_in_ref!(ss, soln.zero)

    populate_solution!(ss)

    if flow_direction == false
        @info "calculated flow direction is opposite to assumed direction in some pipe(s)"
    end

    if length(negative_flow_in_compressors) > 0
        @warn "calculated flow direction is opposite to given direction in some compressor(s)"
        return SolverReturn(compressor_flow_negative, 
            iters_initial + soln.iterations, 
            soln.residual_norm, 
            t_first + t_second, 
            soln.zero, negative_flow_in_compressors)
    end

    return SolverReturn(successfull, 
        iters_initial + soln.iterations, 
        soln.residual_norm, 
        t_first + t_second, 
        soln.zero, Int[])
end

function _create_initial_guess_dof!(ss::SteadySimulator)::Array
    ndofs = length(ref(ss, :dof))
    x_guess = 0.5 * ones(Float64, ndofs) 
    dofs_updated = 0

    components = [:node, :pipe, :compressor, 
        :control_valve, :valve, 
        :resistor, :loss_resistor, :short_pipe]

    for component in components 
        for (i, val) in get(ss.initial_guess, component, [])
            x_guess[ref(ss, component, i, :dof)] = val 
            dofs_updated += 1
        end 
    end 
    return x_guess
end

function check_for_negative_pressures(ss::SteadySimulator, soln_vec::Array)::Bool
    for (i, _) in ref(ss, :node)
        if  soln_vec[ref(ss, :node, i, :dof)] < 0
            return false
        end
    end
    return true
end

function reinitialize_for_positive_pressure!(ss::SteadySimulator, soln_vec::Array)
    for (i, _) in ref(ss, :node)
        val = soln_vec[ref(ss, :node, i, :dof)]
        if  val < 0
            soln_vec[ref(ss, :node, i, :dof)] = abs(val)
        end
    end
end