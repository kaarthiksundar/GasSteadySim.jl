using GasSteadySim
using JSON

file = "./data/GasLib-40-multiple-slacks/"

eos_var = :ideal
ss = initialize_simulator(file, eos=eos_var, initial_guess_filename="")
solver_return = run_simulator!(ss)

println(solver_return.status)

#============== Save solution data for use =============================#
#=
filename = "../test/data/GasLib-40-multiple-slacks/exact_sol_" * string(eos_var) * ".json"
open(filename, "w") do f 
        JSON.print(f, ss.sol, 2)
end
=#
#=======================================================================#
