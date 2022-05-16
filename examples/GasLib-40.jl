using GasSteadySim

file = "./data/GasLib-40/"
ss = initialize_simulator(file, eos=:ideal, initial_guess_filename="")
p = run_simulator!(ss)

# println(solver_return.status)
