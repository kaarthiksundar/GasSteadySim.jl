using GasSteadySim

file = "./data/GasLib-40/"
ss = initialize_simulator(file, eos=:ideal, initial_guess_filename="")
solver_return = run_simulator!(ss, continuation_param=1.0, iteration_limit=20, show_trace=true)


println(solver_return)
