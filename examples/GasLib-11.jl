using GasSteadySim

file = "./data/GasLib-11/"
ss = initialize_simulator(file, eos=:full_cnga, initial_guess_filename="")
solver_return = run_simulator!(ss)

println(solver_return.status)
