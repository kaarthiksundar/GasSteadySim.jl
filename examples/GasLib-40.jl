using GasSteadySim
using JSON

file = "./data/GasLib-40/"
eos_var = :ideal
ss = initialize_simulator(file, eos=eos_var, initial_guess_filename="")
solver_return = run_simulator!(ss, gravity_bool= false, inertial_bool=false, iteration_limit=100, show_trace_flag=true)

println(solver_return.status)


#============== Save solution data for use =============================#
if solver_return.status != nl_solve_failure
        r1 = pipe_equations_no_gravity_no_inertia(ss, solver_return.solution)
        @show r1
        r2 = pipe_equations_no_gravity_with_inertia(ss, solver_return.solution)
        @show r2
        r3 = pipe_equations_with_gravity_no_inertia(ss, solver_return.solution)
        @show r3
        r4 = pipe_equations_with_gravity_with_inertia(ss, solver_return.solution)
        @show r4

        # filename = "./data/GasLib-40/exact_sol_collocation"  * ".json"
        # open(filename, "w") do f 
        #         JSON.print(f, ss.sol, 2)
        # end
end
#======================================================================#