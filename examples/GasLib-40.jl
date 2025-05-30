using GasSteadySim
using JSON

file = "./data/GasLib-40/"
eos_var = :ideal
inertial_bool = false
gravity_bool = false
ss = initialize_simulator(file, eos=eos_var, initial_guess_filename="r.json")
solver_return = run_simulator!(ss, gravity_bool= gravity_bool, inertial_bool=inertial_bool, iteration_limit=100, show_trace_flag=true)

println(solver_return.status)


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

        # filename = "./data/GasLib-40/exact_sol_ideal"  * ".json"
        # open(filename, "w") do f 
        #         JSON.print(f, ss.sol, 2)
        # end
        # filename = "./data/GasLib-40/r4"  * ".json"
        # open(filename, "w") do f 
        #         JSON.print(f, ss.sol, 2)
        # end
end
#======================================================================#