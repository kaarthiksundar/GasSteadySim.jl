using GasSteadySim
using JSON
file = "./data/8-node/"
eos_var = :ideal
ss = initialize_simulator(file, eos=eos_var, initial_guess_filename="")
inertial_bool = false
gravity_bool = false

# ss = initialize_simulator(file, eos=eos_var, initial_guess_filename="exact_sol_ideal.json")

# using IterativeSolvers
# solver_return = run_simulator!(ss, linsolve=(x, A, b) -> copyto!(x, gmres(A, b,  maxiter = 7000) ) )
solver_return = run_simulator!(ss, gravity_bool=gravity_bool, inertial_bool= inertial_bool, method = :newton, iteration_limit = 100, show_trace_flag=true)


println(solver_return.status)
println(ref(ss, :node, 1, "withdrawal") * nominal_values(ss, :mass_flow))


#============== Save solution data for use =============================#
if solver_return.status != nl_solve_failure
        r1 = pipe_equations_no_gravity_no_inertia(ss, solver_return.solution)
        @show r1
        r2 = pipe_equations_no_gravity_with_inertia(ss, solver_return.solution)
        @show r2
        r3 = pipe_equations_with_gravity_no_inertia(ss, solver_return.solution)
        @show r3
        if gravity_bool == true
                r4 = pipe_equations_with_gravity_with_inertia(ss, solver_return.solution)
                @show r4
        end
        # @show r1, r3
        # @show r1, r2, r3, r4
        # filename = "data/8-node/exact_sol_" * string(eos_var) * ".json"
        # open(filename, "w") do f 
        #         JSON.print(f, ss.sol, 2)
        # end
        # @info "Saved result file $filename"
end

#=======================================================================#




 

