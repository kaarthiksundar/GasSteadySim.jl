# helper functions to separate out prep and solve 
function prepare_for_solve!(ss::SteadySimulator, case::Symbol)::NonlinearFunction
    residual_fun! = (r_dof, x_dof, _) -> assemble_residual!(ss, x_dof, r_dof, case) # added _
    Jacobian_fun! = (J_dof, x_dof, _) -> assemble_mat!(ss, x_dof, J_dof, case) # added _
    n = length(ref(ss, :dof))
    J0 = spzeros(n, n)
    assemble_mat!(ss, rand(n), J0, case)
    df = NonlinearFunction(residual_fun!; jac = Jacobian_fun!, jac_prototype = J0)
    return df
end


# run_simulator
function run_simulator!(ss::SteadySimulator; 
    gravity_bool::Bool=false, 
    inertial_bool::Bool=false,
    method::Symbol=:trust_region,
    iteration_limit::Int64=2000, 
    show_trace_flag::Bool=false,
    x_guess::Vector=Vector{Float64}(), 
    kwargs...)::SolverReturn
    
    ss.params[:inertial_bool] = inertial_bool
    ss.params[:gravity_bool] = gravity_bool
    fcn_method = get(solver_method, method, TrustRegion())

    (isempty(x_guess)) && (x_guess = _create_initial_guess_dof!(ss))
    x_guess_original = copy(x_guess)

    df_tpc = prepare_for_solve!(ss, :two_point_collocation)
    prob = NonlinearProblem(df_tpc, x_guess)
    @info "Solving system after formulation of pipe equation as Two-Point Collocation..."
    time = @elapsed soln = solve(prob, fcn_method; maxiters = iteration_limit, 
    show_trace = Val(show_trace_flag),kwargs...)
    res = maximum(abs.(soln.resid))

    if ss.params[:eos] == :ideal
        err_max, err_rms = get_applicable_residual(ss, soln.u)
        @info "Two-Point Collocation residual for ideal EoS integral - \n max: $err_max, rms: $err_rms"
    end


    if  res < 1e-4 || SciMLBase.successful_retcode(soln) == false
        @info "Using solution from Two-Point Collocation as initial guess..."
        x_guess = soln.u
     else
        @info "NOT using solution from Two-Point Collocation. Using given initial guess..."
        x_guess = x_guess_original
    end

    df_ode = prepare_for_solve!(ss, :ode)
    prob = NonlinearProblem(df_ode, x_guess)
    @info "Solving system after formulation of pipe equation as an ODE..."

    time = @elapsed soln = solve(prob, fcn_method; maxiters = iteration_limit, 
    show_trace = Val(show_trace_flag),kwargs...)
    res = maximum(abs.(soln.resid))

    if ss.params[:eos] == :ideal
        err_max, err_rms = get_applicable_residual(ss, soln.u)
        @info "ODE formulation residual for ideal EoS integral - \n max: $err_max, rms: $err_rms"
    end


    convergence_state = SciMLBase.successful_retcode(soln) || (res < 1e-4) # do this correctly
    @info soln.retcode
    
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

function _create_initial_guess_dof!(ss::SteadySimulator)::Array
    ndofs = length(ref(ss, :dof))
    Random.seed!(2025)
    x_guess = rand(ndofs) 

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