using GasSteadySim
using JSON

file = "./data/GasLib-40/"
eos_var = :ideal
inertial_bool = false
gravity_bool = false
ss = initialize_simulator(file, eos=eos_var, initial_guess_filename="r1.json")
solver_return = run_simulator!(ss, method=:trust_region, gravity_bool= gravity_bool, inertial_bool=inertial_bool, iteration_limit=100, show_trace_flag=false, reltol = 1e-3)
println(solver_return.status)


#============== Save solution data for use =============================#
if solver_return.status != nl_solve_failure
        # filename = "./data/GasLib-40/exact_sol_ideal"  * ".json"
        # open(filename, "w") do f 
        #         JSON.print(f, ss.sol, 2)
        # end
end
#======================================================================#