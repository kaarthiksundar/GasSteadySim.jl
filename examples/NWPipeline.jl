using GasSteadySim
import JSON

file = "./data/NWPipeline/"
eos_var = :simple_cnga
ss = initialize_simulator(file, eos=eos_var, initial_guess_filename="NR.json")
# ss = initialize_simulator(file, eos=eos_var, case_name= "one", case_types =[:bc, :network],initial_guess_filename="NR.json")

solver_return = run_simulator!(ss; method=:trust_region, show_trace_flag=true)

println(solver_return.status)

#============== Save solution data for use =============================#

# filename = file * "exact_sol_$eos_var.json"
filename = file * "NR.json"
open(filename, "w") do f
        JSON.print(f, ss.sol, 2)
end

#=======================================================================#
