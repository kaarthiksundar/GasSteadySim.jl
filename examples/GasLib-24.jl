using GasSteadyODE
using JSON

file = "./data/GasLib-24/"
eos_var = :ideal
inertial_bool = false
gravity_bool = false
guess_file = "exact_sol_ideal.json"
write_bool = false
write_file = ""

ss = initialize_simulator(file, eos=eos_var, initial_guess_filename=guess_file)
solver_return = run_simulator!(ss, method=:newton, gravity_bool= gravity_bool, inertial_bool=inertial_bool, iteration_limit=100, show_trace_flag=true, reltol = 1e-3)

println(solver_return.status)


#============== Save solution data for use =============================#
if solver_return.status != nl_solve_failure
        if write_bool
            filename = file * write_file
            open(filename, "w") do f 
                JSON.print(f, ss.sol, 2)
            end
        end
        
end
#======================================================================#