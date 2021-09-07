using GasSteadySim

file = "./data/GasLib-40/"
ss = initialize_simulator(file, eos=:simple_cnga, initial_guess_filename="")
solver_return = run_simulator!(ss)

println(solver_return.status)
