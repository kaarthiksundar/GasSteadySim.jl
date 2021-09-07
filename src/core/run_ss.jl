function run_simulator!(ss::SteadySimulator; 
    jacobian_type::Symbol=:in_place, 
    method::Symbol=:trust_region,
    iteration_limit::Int64=2000)::SolverReturn
    
    x_guess = _create_initial_guess_dof!(ss)
    
    if jacobian_type == :sparse
        residual_fun = (x_dof) -> assemble_residual(ss, x_dof)
        Jacobian_fun = (x_dof) -> assemble_mat(ss, x_dof)
        J0 = Jacobian_fun(x_guess)
        df = OnceDifferentiable(residual_fun, Jacobian_fun, x_guess, x_guess, J0)
    end

    if jacobian_type == :in_place
        residual_fun = (r_dof, x_dof) -> assemble_residual_in_place(ss, x_dof, r_dof)
        Jacobian_fun = (J_dof, x_dof) -> assemble_mat_in_place(ss, x_dof, J_dof)
        df = nothing
    end

    t_first, soln = _invoke_solver(residual_fun, Jacobian_fun, df, x_guess, method, iteration_limit)

    t_second = 0
    all_pressures_non_neg = check_for_negative_pressures(ss, soln.zero)
    convergence_state = converged(soln)

    if convergence_state == false
        return SolverReturn(initial_nl_solve_failure, 
            soln.iterations, 
            soln.residual_norm, 
            t_first + t_second, 
            soln.zero)
    end

    if all_pressures_non_neg == true
        @info "initial nonlinear solve complete and all pressures non-negative"
    end

    iters_initial = soln.iterations

    if all_pressures_non_neg == false
        @info "correcting pressures..."
        reinitialize_for_positive_pressure!(ss, soln.zero)
        t_second, soln = _invoke_solver(residual_fun, Jacobian_fun, df, soln.zero, method, iteration_limit)
        
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
                soln.zero)
        end

        if all_pressures_non_neg == false
            @warn "pressure correction failed"
            update_solution_fields_in_ref!(ss, soln.zero)
            populate_solution!(ss)
            return SolverReturn(pressure_correction_failure, 
                iters_initial + soln.iterations, 
                soln.residual_norm, 
                t_first + t_second, 
                soln.zero)
        else 
            @info "pressure correction successful"
        end
    end


    flow_direction, compressor_direction = update_solution_fields_in_ref!(ss, soln.zero)

    populate_solution!(ss)

    if flow_direction == false
        @info "calculated flow direction is opposite to assumed direction in some pipe(s)"
    end

    if compressor_direction == false
        @warn "calculated flow direction is opposite to given direction in some compressor(s)"
        return SolverReturn(compressor_flow_negative, 
            iters_initial + soln.iterations, 
            soln.residual_norm, 
            t_first + t_second, 
            soln.zero)
    end

    return SolverReturn(successfull, 
        iters_initial + soln.iterations, 
        soln.residual_norm, 
        t_first + t_second, 
        soln.zero)
end

function _invoke_solver(r, J, df, x_0, method, iter)
    if isa(df, Nothing)
        t = @elapsed soln = nlsolve(r, J, x_0; method = method, iterations = iter)
    else 
        t = @elapsed soln = nlsolve(df, x_0; method = method, iterations = iter)
    end 
    return t, soln
end 

function _create_initial_guess_dof!(ss::SteadySimulator)::Array
    ndofs = length(ref(ss, :dof))
    x_guess = 0.5 * ones(Float64, ndofs) 
    dofs_updated = 0

    for (i, val) in get(ss.initial_guess, :node, [])
        x_guess[ref(ss, :node, i, :dof)] = val
        dofs_updated += 1
    end

    for (i, val) in get(ss.initial_guess, :pipe, [])
        x_guess[ref(ss, :pipe, i, :dof)] = val
        dofs_updated += 1 
    end

    for (i, val) in get(ss.initial_guess, :compressor, [])
        x_guess[ref(ss, :compressor, i, :dof)] = val
        dofs_updated += 1
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