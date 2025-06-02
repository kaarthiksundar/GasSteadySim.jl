# helper functions to separate out prep and solve 
function prepare_for_solve!(ss::SteadySimulator)::NonlinearFunction
    residual_fun! = (r_dof, x_dof, _) -> assemble_residual!(ss, x_dof, r_dof) # added _
    Jacobian_fun! = (J_dof, x_dof, _) -> assemble_mat!(ss, x_dof, J_dof) # added _
    n = length(ref(ss, :dof))
    J0 = spzeros(n, n)
    assemble_mat!(ss, rand(n), J0)
    df = NonlinearFunction(residual_fun!; jac = Jacobian_fun!, jac_prototype = J0)
    # df = OnceDifferentiable(residual_fun!, Jacobian_fun!, rand(n), rand(n), J0)
    return df
end

# run function
function run_simulator!(
    ss::SteadySimulator, 
    df::NonlinearFunction; 
    x_guess::Vector=Vector{Float64}(), 
    method::Symbol=:trust_region,
    iteration_limit::Int64=2000, 
    show_trace_flag::Bool=false,
    kwargs...)::SolverReturn

    (isempty(x_guess)) && (x_guess = _create_initial_guess_dof!(ss))
    prob = NonlinearProblem(df, x_guess)
    fcn_method = get(solver_method, method, TrustRegion())
    time = @elapsed soln = solve(prob, fcn_method, 
    maxiters = iteration_limit, show_trace = Val(show_trace_flag))
    res = maximum(abs.(soln.resid), kwargs...)
    convergence_state = soln.retcode

    if convergence_state == false
        return SolverReturn(nl_solve_failure, 
            soln.stats, 
            res, 
            time, soln.u, 
            Int[], Int[])
    end

    sol_return = update_solution_fields_in_ref!(ss, soln.u)
    populate_solution!(ss)

    unphysical_solution_flag = ~isempty(sol_return[:compressors_with_neg_flow]) || ~isempty(sol_return[:nodes_with_negative_pressures])

    if unphysical_solution_flag

        return SolverReturn(unphysical_solution, 
                soln.stats, 
                res, 
                time, soln.u, 
                sol_return[:compressors_with_neg_flow], 
                sol_return[:nodes_with_negative_pressures])
        
    end 

    return SolverReturn(physical_solution, 
        soln.stats, 
        res, 
        time, soln.u, 
        sol_return[:compressors_with_neg_flow], 
        sol_return[:nodes_with_negative_pressures])
end

# overloaded run_simulator
function run_simulator!(ss::SteadySimulator; 
    gravity_bool::Bool=false, 
    inertial_bool::Bool=false,
    method::Symbol=:trust_region,
    iteration_limit::Int64=2000, 
    show_trace_flag::Bool=false,
    kwargs...)::SolverReturn

    ss.params[:inertial_bool] = inertial_bool
    ss.params[:gravity_bool] = gravity_bool
    
    x_guess = _create_initial_guess_dof!(ss)
    df = prepare_for_solve!(ss)
    return run_simulator!(ss, df, 
        x_guess = x_guess, 
        method = method, 
        iteration_limit = iteration_limit, 
        show_trace_flag = show_trace_flag, kwargs...)
end

function _create_initial_guess_dof!(ss::SteadySimulator)::Array
    ndofs = length(ref(ss, :dof))
    Random.seed!(2025)
    x_guess = rand(ndofs)
    # x_guess = 0.5 * ones(Float64, ndofs) 
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