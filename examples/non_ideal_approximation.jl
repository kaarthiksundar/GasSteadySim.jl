using GasSteadySim
using JSON
using LinearAlgebra
using NLSolversBase

file = "./data/GasLib-11/"


function create_cnga_solution_pressure_formulation(
    ss_pressure::SteadySimulator, 
    ss_potential::SteadySimulator)::Vector{Float64}
    
    ndofs = length(ref(ss_pressure, :dof))
    x_guess = zeros(Float64, ndofs) 

    for i = 1:length(ref(ss_potential, :dof))
        comp, id = ss_potential.ref[:dof][i]
        if comp == :node
            val = ref(ss_pressure, :is_pressure_node, id) ? ss_potential.ref[comp][id]["pressure"] : ss_potential.ref[comp][id]["potential"]
            x_guess[ss_pressure.ref[comp][id]["dof"]] = val
        elseif comp == :pipe
            x_guess[ss_pressure.ref[comp][id]["dof"]] = ss_potential.ref[comp][id]["flow"]
        elseif comp == :compressor
            x_guess[ss_pressure.ref[comp][id]["dof"]] = ss_potential.ref[comp][id]["flow"]
        elseif comp == :valve 
            x_guess[ss_pressure.ref[comp][id]["dof"]] = ss_potential.ref[comp][id]["flow"]
        end
    end

    return x_guess
end


eos = :simple_cnga
ss_potential = initialize_simulator(file, 
    eos=eos, 
    use_potential_formulation=true, 
    potential_ratio_coefficients=[0.0, 0.0, 0.9, 0.1], 
    initial_guess_filename="") 

df_potential = prepare_for_solve!(ss_potential)
solver_potential = run_simulator!(ss_potential, df_potential)

ss_pressure = initialize_simulator(file, 
    eos=eos, 
    initial_guess_filename="") 

df_pressure = prepare_for_solve!(ss_pressure)

x_guess = create_cnga_solution_pressure_formulation(ss_pressure, ss_potential)

num_compressors = length(ss_pressure.ref[:compressor])
var = value!(df_pressure, x_guess)

@show collect(var)

println("Pipe + Node Residual (inf norm): ", norm(var[1:end-num_compressors], Inf), "\nCompressor  Residual (inf norm): ", norm(var[end-num_compressors+1:end], Inf))

# println(findall(x->abs(x)>1e-3, var)," ", var, " ", norm(var), " ", norm(var, Inf))

solver_pressure = run_simulator!(ss_pressure, df_pressure, x_guess = x_guess, show_trace_flag=true)
# solver_p = solve_on_network!(ss, df, x_guess=x_dof, iteration_limit=1)

println(solver_pressure.iterations, " ", solver_pressure.residual_norm)
println("Absolute error of approx soln (inf norm): ", norm(x_guess - solver_pressure.solution, Inf), "\nRelative error of approximate soln (inf norm): ", norm(x_guess - solver_pressure.solution, Inf)/norm(solver_pressure.solution, Inf))



# for i = 1:length(ss_p.ref[:node])
#     v1 = ss_p.ref[:node][i]["potential"]
#     v2 = ss_pot.ref[:node][i]["potential"]
#     println(v1, " ", v2, " ", abs(v1-v2))
# end
# println("nodes finished")
# for i = 1:length(ss_p.ref[:pipe])
#     v1 = ss_p.ref[:pipe][i]["flow"]
#     v2 = ss_pot.ref[:pipe][i]["flow"]
#     println(v1, " ", v2, " ", abs(v1-v2))
# end
# println("pipes finished")

# for i = 1:length(ss_p.ref[:compressor])
#     v1 = ss_p.ref[:compressor][i]["flow"]
#     v2 = ss_pot.ref[:compressor][i]["flow"]
#     println(v1, " ", v2, " ", abs(v1-v2))
# end
# println("compressors finished")

