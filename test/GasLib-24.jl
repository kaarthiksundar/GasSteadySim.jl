

@testset "test GasLib-24 ideal run" begin
    file = "./data/GasLib-24/"
    ss = initialize_simulator(file, eos=:ideal, initial_guess_filename="exact_sol_ideal.json")
    solver_return = run_simulator!(ss)

    @show solver_return.stats
    exact_sol = GasSteadyODE._parse_json(file * "exact_sol_ideal.json")
    _check_correctness(ss.sol, exact_sol)


end

@testset "test GasLib-24 simple CNGA run" begin
    file = "./data/GasLib-24/"
    ss = initialize_simulator(file, eos=:simple_cnga, initial_guess_filename="exact_sol_simple_cnga")
    solver_return = run_simulator!(ss)

    @show solver_return.stats
    exact_sol = GasSteadyODE._parse_json(file * "exact_sol_simple_cnga.json")
    _check_correctness(ss.sol, exact_sol)
end


@testset "test GasLib-24 full CNGA run" begin
    file = "./data/GasLib-24/"
    ss = initialize_simulator(file, eos=:full_cnga, initial_guess_filename="exact_sol_full_cnga")
    solver_return = run_simulator!(ss)

    @show solver_return.stats
    exact_sol = GasSteadyODE._parse_json(file *"exact_sol_full_cnga.json")
    _check_correctness(ss.sol, exact_sol)
end
