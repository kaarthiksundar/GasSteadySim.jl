using GasSteadyODE
using JSON
file = "./data/8-node/"
eos_var = :ideal

# using IterativeSolvers
# solver_return = run_simulator!(ss, linsolve=(x, A, b) -> copyto!(x, gmres(A, b,  maxiter = 7000) ) )
inertial_bool = false
gravity_bool = false
guess_file = "r1.json"
write_bool = false
write_file = "r4.json"

ss = initialize_simulator(file, eos=eos_var, initial_guess_filename=guess_file)
solver_return = run_simulator!(ss, method=:trust_region, gravity_bool= gravity_bool, inertial_bool=inertial_bool, iteration_limit=100, show_trace_flag=true, reltol = 1e-3)

println(solver_return.status)


#============== Save solution data for use =============================#
if solver_return.status != nl_solve_failure
        if write_bool
            filename = "./data/8-node/" * write_file
            open(filename, "w") do f 
                JSON.print(f, ss.sol, 2)
            end
        end
        
end
#===========================#

