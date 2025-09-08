# helper functions to separate out prep and solve 
function prepare_for_solve!(ss::SteadySimulator)::NonlinearFunction
    residual_fun! = (r_dof, x_dof, _) -> assemble_residual!(ss, x_dof, r_dof)
    Jacobian_fun! = (J_dof, x_dof, _) -> assemble_mat!(ss, x_dof, J_dof)
    n = length(ref(ss, :dof))
    J0 = spzeros(n, n)
    assemble_mat!(ss, rand(n), J0)
    df = NonlinearFunction(residual_fun!; jac = Jacobian_fun!, jac_prototype = J0)
    return df
end

# run function
function run_simulator!(
    ss::SteadySimulator;
    x_guess::Vector=Vector{Float64}(), 
    method::Symbol=:trust_region,
    iteration_limit::Int64=2000, 
    show_trace_flag::Bool=false,
    kwargs...)::SolverReturn

    (isempty(x_guess)) && (x_guess = _create_initial_guess_dof!(ss))
    fcn_method = get(solver_method, method, TrustRegion())
    df = prepare_for_solve!(ss)
    prob = NonlinearProblem(df, x_guess)
    time = @elapsed soln = solve(prob, fcn_method; maxiters = iteration_limit, 
        show_trace = Val(show_trace_flag), kwargs...)
    res = maximum(abs.(soln.resid))
    
    convergence_state = SciMLBase.successful_retcode(soln) == true

    if convergence_state == false
        return SolverReturn(nl_solve_failure, 
            soln.stats, 
            res, 
            time, soln.u, 
            Int[], Int[], Int[])
    end

    sol_return = update_solution_fields_in_ref!(ss, soln.u)
    populate_solution!(ss)

    is_solution_unphysical = ~isempty(sol_return[:compressors_with_neg_flow]) || 
    ~isempty(sol_return[:nodes_with_neg_potential])

    if is_solution_unphysical
        is_unique = isempty(sol_return[:nodes_with_pressure_not_in_domain])
        if is_unique 
            return SolverReturn(unique_unphysical_solution, 
                soln.stats, 
                res, 
                time, soln.u, 
                sol_return[:compressors_with_neg_flow], 
                sol_return[:nodes_with_neg_potential],
                sol_return[:nodes_with_pressure_not_in_domain])
        else 
            return SolverReturn(unphysical_solution, 
                soln.stats, 
                res, 
                time, soln.u, 
                sol_return[:compressors_with_neg_flow], 
                sol_return[:nodes_with_neg_potential],
                sol_return[:nodes_with_pressure_not_in_domain])
        end 
    end 

    return SolverReturn(unique_physical_solution, 
        soln.stats, 
        res, 
        time, soln.u, 
        sol_return[:compressors_with_neg_flow], 
        sol_return[:nodes_with_neg_potential],
        sol_return[:nodes_with_pressure_not_in_domain])
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
            if isa(val,Nothing) 
                continue
            end
            x_guess[ref(ss, component, i, "dof")] = val 
            dofs_updated += 1
        end 
    end
    if 0 < dofs_updated < ndofs
        @warn "Null values found in ig file replaced by 0.5"
    end 
    return x_guess
end