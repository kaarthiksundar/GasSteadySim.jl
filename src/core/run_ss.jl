function run_simulator!(ss::SteadySimulator; 
    method::Symbol=:newton,
    iteration_limit::Int64=2000, 
    kwargs...)::SolverReturn
    
    x_guess = _create_initial_guess_dof!(ss)
    n = length(x_guess)
    

    residual_fun! = (r_dof, x_dof) -> assemble_residual!(ss, x_dof, r_dof)
    Jacobian_fun! = (J_dof, x_dof) -> assemble_mat!(ss, x_dof, J_dof)
    J0 = spzeros(n, n)
    assemble_mat!(ss, rand(n), J0)
    df = OnceDifferentiable(residual_fun!, Jacobian_fun!, rand(n), rand(n), J0)

    t_first = @elapsed soln = nlsolve(df, x_guess; method = method, iterations = iteration_limit, kwargs...)

    t_second = 0.0
    convergence_state = converged(soln)
    pressure_correction_performed = false 
    compressor_pressures_in_proper_range = check_compressor_pressures(ss, soln.zero)
    neg_potential_exists = check_for_negative_potentials(ss, soln.zero)
    all_pressures_non_neg = check_for_negative_pressures(ss, soln.zero)

    if convergence_state == false
        return SolverReturn(initial_nl_solve_failure, 
            pressure_correction_performed,
            soln.iterations, 
            soln.residual_norm, 
            t_first + t_second, 
            soln.zero, Int[])
    end

    _, negative_flow_in_compressors = update_solution_fields_in_ref!(ss, soln.zero)

    if compressor_pressures_in_proper_range == false 
        populate_solution!(ss)
        return SolverReturn(compressor_pressure_not_in_proper_range,
            pressure_correction_performed,
            soln.iterations, 
            soln.residual_norm, 
            t_first + t_second, 
            soln.zero, Int[])
    end 

    if length(negative_flow_in_compressors) > 0
        @warn "calculated flow direction is opposite to given direction in some compressor(s)"
        populate_solution!(ss)
        return SolverReturn(compressor_flow_negative, 
            pressure_correction_performed,
            soln.iterations, 
            soln.residual_norm, 
            t_first + t_second, 
            soln.zero, Int[])
    end

    if neg_potential_exists == true 
        populate_solution!(ss)
        return SolverReturn(potential_negative, 
        pressure_correction_performed,
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
        pressure_correction_performed = true 
        reinitialize_for_positive_pressure!(ss, soln.zero)
        t_second = @elapsed soln = nlsolve(df, soln.zero; method = method, iterations = iteration_limit)
        
        all_pressures_non_neg = check_for_negative_pressures(ss, soln.zero)
        convergence_state = converged(soln)

        if convergence_state == false
            @warn "nonlinear solve for pressure correction failed!"
            update_solution_fields_in_ref!(ss, soln.zero)
            populate_solution!(ss)
            return SolverReturn(pressure_correction_nl_solve_failure, 
                pressure_correction_performed,
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
                pressure_correction_performed,
                iters_initial + soln.iterations, 
                soln.residual_norm, 
                t_first + t_second, 
                soln.zero, Int[])
        else 
            @info "pressure correction successful"
        end
    end

    return SolverReturn(successfull, 
        pressure_correction_performed,
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
            x_guess[ref(ss, component, i, "dof")] = val 
            dofs_updated += 1
        end 
    end 
    return x_guess
end

function check_compressor_pressures(ss::SteadySimulator, soln_vec::Array)::Bool 
    b1, b2 = get_eos_coeffs(ss) 
    lb = (b2 == 0.0) ? 0.0 : (- 3.0) * b1 / 2.0 / b2
    for (_, compressor) in ref(ss, :compressor)
        fr_node = compressor["fr_node"]
        to_node = compressor["to_node"]
        fr_p = soln_vec[ref(ss, :node, fr_node, "dof")]
        to_p = soln_vec[ref(ss, :node, to_node, "dof")]
        fr = (lb == 0.0) ? (fr_p < 0.0) : (fr_p > lb && fr_p < 0.0)
        to = (lb == 0.0) ? (to_p < 0.0) : (to_p > lb && to_p < 0.0)
        (fr || to) && (return false)
    end 
    return true
end 

function check_for_negative_potentials(ss::SteadySimulator, soln_vec::Array)::Bool 
    b1, b2 = get_eos_coeffs(ss) 
    for (i, _) in ref(ss, :node)
        p = soln_vec[ref(ss, :node, i, "dof")]
        potential = b1 / 2.0 * p^2 + b2 / 3.0 * p^3
        (potential < 0.0) && (return true)
    end 
    return false
end 

function check_for_negative_pressures(ss::SteadySimulator, soln_vec::Array)::Bool
    for (i, _) in ref(ss, :node)
        if  soln_vec[ref(ss, :node, i, "dof")] < 0
            return false
        end
    end
    return true
end

function reinitialize_for_positive_pressure!(ss::SteadySimulator, soln_vec::Array)
    for (i, _) in ref(ss, :node)
        val = soln_vec[ref(ss, :node, i, "dof")]
        if  val < 0
            soln_vec[ref(ss, :node, i, "dof")] = abs(val)
        end
    end
end