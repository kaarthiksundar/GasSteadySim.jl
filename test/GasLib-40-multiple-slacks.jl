

@testset "test GasLib-40-multiple-slacks ideal run" begin
    file = "./data/GasLib-40-multiple-slacks/"
    ss = initialize_simulator(file, eos=:ideal, initial_guess_filename="")
    solver_return = run_simulator!(ss)

    @show solver_return.iterations
    exact_sol = GasSteadySim._parse_json("data/GasLib-40-multiple-slacks/exact_sol_ideal.json")
    _check_correctness(ss.sol, exact_sol)


end

@testset "test GasLib-40-multiple-slacks simple CNGA run" begin
    file = "./data/GasLib-40-multiple-slacks/"
    ss = initialize_simulator(file, eos=:simple_cnga, initial_guess_filename="")
    solver_return = run_simulator!(ss)

    @show solver_return.iterations
    exact_sol = GasSteadySim._parse_json("data/GasLib-40-multiple-slacks/exact_sol_simple_cnga.json")
    _check_correctness(ss.sol, exact_sol)
end


@testset "test GasLib-40-multiple-slacks full CNGA run" begin
    file = "./data/GasLib-40-multiple-slacks/"
    ss = initialize_simulator(file, eos=:full_cnga, initial_guess_filename="")
    solver_return = run_simulator!(ss)

    @show solver_return.iterations
    exact_sol = GasSteadySim._parse_json("data/GasLib-40-multiple-slacks/exact_sol_full_cnga.json")
    _check_correctness(ss.sol, exact_sol)
end
