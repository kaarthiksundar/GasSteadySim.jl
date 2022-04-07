function run_simulator!(ss::SteadySimulator; 
    method::Symbol=:newton,
    iteration_limit::Int64=2000, 
    kwargs...)::SolverReturn
    
    x_guess = _create_initial_guess_dof!(ss)
    n = length(x_guess)

    # contraction_factor = 1/100
    # r_dof =  zeros(Float64, n) 
    # residual_fun_fixed_point! = (r_dof, x_dof) -> assemble_residual_fixed_point_iteration!(ss, x_dof, r_dof, contraction_factor)
    # for iter = 1:100
    #     residual_fun_fixed_point!(r_dof, x_guess)
    #     # println(x_guess, "\n", r_dof, "\n", norm(r_dof - x_guess), "\n")
    #     println(norm(r_dof - x_guess), "\n")

    #     for i =1:n
    #         x_guess[i] = r_dof[i]
    #     end
    # end
    # println(x_guess, "\n")
    # # return SolverReturn(successfull, 4, norm(r_dof), 0.001, x_guess, Int[])
    

    residual_fun! = (r_dof, x_dof) -> assemble_residual!(ss, x_dof, r_dof)
    Jacobian_fun! = (J_dof, x_dof) -> assemble_mat!(ss, x_dof, J_dof)
    J0 = spzeros(n, n)
    assemble_mat!(ss, rand(n), J0)
    df = OnceDifferentiable(residual_fun!, Jacobian_fun!, rand(n), rand(n), J0)

    t_first = @elapsed soln = nlsolve(df, x_guess; method = method, iterations = iteration_limit, kwargs...)

    convergence_state = converged(soln)


    if convergence_state == false
        return SolverReturn(nl_solve_failure, 
            soln.iterations, 
            soln.residual_norm, 
            t_first, 
            soln.zero, Int[], Int[])
    end

    _, negative_flow_in_compressors, negative_nodal_pressures, nodal_pressures_not_in_domain = update_solution_fields_in_ref!(ss, soln.zero)

    println(negative_flow_in_compressors, negative_nodal_pressures, nodal_pressures_not_in_domain, "\n")
    populate_solution!(ss)


    if ~isempty(negative_flow_in_compressors) || ~isempty(negative_nodal_pressures) == true

        if ~isempty(nodal_pressures_not_in_domain)
            @warn "unphysical solution"
            return SolverReturn(unphysical_soln, 
                soln.iterations, 
                soln.residual_norm, 
                t_first,
                soln.zero, negative_flow_in_compressors, negative_nodal_pressures)
        else
            @warn "unphysical, unique solution"
            return SolverReturn(unique_unphysical_soln, 
                soln.iterations, 
                soln.residual_norm, 
                t_first, 
                soln.zero, negative_flow_in_compressors, negative_nodal_pressures)
        end 
    end 

    
    return SolverReturn(unique_feasible_soln, 
        soln.iterations, 
        soln.residual_norm, 
        t_first, 
        soln.zero, negative_flow_in_compressors, negative_nodal_pressures)
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

# # checks for hypothesis 1 in the paper
# function _check_pressure_hypothesis(ss::SteadySimulator, soln_vec::Array)::Bool 
#     b1, b2 = get_eos_coeffs(ss) 
#     lb = (b2 == 0.0) ? 0.0 : (- 3.0) * b1 / 2.0 / b2
#     for (_, compressor) in ref(ss, :compressor)
#         fr_node = compressor["fr_node"]
#         to_node = compressor["to_node"]
#         fr_p = soln_vec[ref(ss, :node, fr_node, "dof")]
#         to_p = soln_vec[ref(ss, :node, to_node, "dof")]
#         fr = (lb == 0.0) ? (fr_p < 0.0) : (fr_p > lb && fr_p < 0.0)
#         to = (lb == 0.0) ? (to_p < 0.0) : (to_p > lb && to_p < 0.0)
#         (fr || to) && (return false)
#     end 
#     return true
# end 

# function _check_for_negative_potentials_and_pressures(ss::SteadySimulator, soln_vec::Array)::Tuple{Bool, Bool} 
    
#     physical_flag = true
#     uniqueness_flag = true

#     for (i, _) in ref(ss, :node)
#         if ref(ss, :is_pressure_node, i)
#             pressure = soln_vec[ref(ss, :node, i, "dof")
#             potential = get_potential(ss, pressure)
#             if pressure < 0.0 && potential 
#         else 
#             potential = soln_vec[ref(ss, :node, i, "dof")]
#         end
#         (potential < 0.0) && (return true)
#     end 
#     return false
# end 

# function _check_for_negative_pressures(ss::SteadySimulator, soln_vec::Array)::Bool
#     for (i, _) in ref(ss, :node)
#         if  soln_vec[ref(ss, :node, i, "dof")] < 0
#             return false
#         end
#     end
#     return true
# end

# function reinitialize_for_positive_pressure!(ss::SteadySimulator, soln_vec::Array)
#     for (i, _) in ref(ss, :node)
#         val = soln_vec[ref(ss, :node, i, "dof")]
#         if  val < 0
#             soln_vec[ref(ss, :node, i, "dof")] = abs(val)
#         end
#     end
# end