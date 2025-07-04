using GasSteadySim
import JSON

file = "./data/GasLib-40/"
eos_var = :full_cnga
ss = initialize_simulator(file, eos=eos_var, initial_guess_filename="")
solver_return = run_simulator!(ss)

println(solver_return.status)

#============== Save solution data for use =============================#

# filename = file * "exact_sol_$eos_var.json"
# open(filename, "w") do f 
#         JSON.print(f, ss.sol, 2)
# end

#=======================================================================#
