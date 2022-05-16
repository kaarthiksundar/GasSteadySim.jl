using GasSteadySim





file = "./data/8-node/"
ss = initialize_simulator(file, eos=:ideal, initial_guess_filename="")
p = run_simulator!(ss)

# println(solver_return.status)
# println(ref(ss, :node, 1, "withdrawal") * nominal_values(ss, :mass_flow))


#

# x = Variable()
# xmax = 1
# y = Semidefinite(num_ns)
# z = Semidefinite(num_edges)
# # y = Variable(num_ns,  Positive() )
# # z = Variable(num_edges,  Positive() )

# # SDP constraints
# # p = minimize(-x, diagm(y + x) * A_ns - A_alpha_ns * diagm(z + x)    == zeros(num_ns, num_edges), x > 0,  x <= xmax )
# p = minimize(-x, (y + x) * A_ns - A_alpha_ns * (z + x)   == zeros(num_ns, num_edges), x > 0,  x <= xmax, opnorm(y, 2) <= xmax, opnorm(z, 2) <= xmax )
# solve!(p, SCS.Optimizer; silent_solver = false)
# println(p.status)


# # println(evaluate(y), "\n", evaluate(z))
# println( det( evaluate( y + x * I(num_ns) )), "\n", det( evaluate( z + x * I(num_edges) ) ) )




# # Generate random problem data
# m = 4;  n = 5
# A = randn(m, n); b = randn(m, 1)

# # Create a (column vector) variable of size n x 1.
# x = Variable(n)

# # The problem is to minimize ||Ax - b||^2 subject to x >= 0
# # This can be done by: minimize(objective, constraints)
# problem = minimize(sumsquares(A * x - b), [x >= 0])

# # Solve the problem by calling solve!
# solve!(problem, SCS.Optimizer; silent_solver = true)

# # Check the status of the problem
# problem.status # :Optimal, :Infeasible, :Unbounded etc.

# # Get the optimum value
# problem.optval

