using GasSteadySim

file = "./data/8-node/"
ss = initialize_simulator(file, eos=:ideal, initial_guess_filename="")
solver_return = run_simulator!(ss)

println(solver_return.status)
println(ref(ss, :node, 1, "withdrawal") * nominal_values(ss, :mass_flow))
