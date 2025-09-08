

@testset "test GasLib-135 ideal run" begin
    file = "./data/GasLib-135/"
    ss = initialize_simulator(file, eos=:ideal, initial_guess_filename="")
    solver_return = run_simulator!(ss)

    @show solver_return.stats.nsteps
    exact_sol = GasSteadySim._parse_json(file * "exact_sol_ideal.json")
    _check_correctness(ss.sol, exact_sol, 5e-4)


end

@testset "test GasLib-135 simple CNGA run" begin
    file = "./data/GasLib-135/"
    ss = initialize_simulator(file, eos=:simple_cnga, initial_guess_filename="")
    solver_return = run_simulator!(ss)

    @show solver_return.stats.nsteps
    exact_sol = GasSteadySim._parse_json(file * "exact_sol_simple_cnga.json")
    _check_correctness(ss.sol, exact_sol, 5e-4)
end


@testset "test GasLib-135 full CNGA run" begin
    file = "./data/GasLib-135/"
    ss = initialize_simulator(file, eos=:full_cnga, initial_guess_filename="")
    solver_return = run_simulator!(ss)

    @show solver_return.stats.nsteps
    exact_sol = GasSteadySim._parse_json(file * "exact_sol_full_cnga.json")
    _check_correctness(ss.sol, exact_sol, 5e-4)
end
