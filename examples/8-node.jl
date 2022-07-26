using GasSteadySim
file = "./data/8-node/"
ss = initialize_simulator(file, eos=:ideal, initial_guess_filename="")

# using IterativeSolvers
# solver_return = run_simulator!(ss, linsolve=(x, A, b) -> copyto!(x, gmres(A, b,  maxiter = 7000) ) )


solver_return = run_simulator!(ss, continuation_param=1.0, iteration_limit=20, show_trace=true)


println(solver_return)
println(ref(ss, :node, 1, "withdrawal") * nominal_values(ss, :mass_flow))
