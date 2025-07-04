@testset "test GasLib-11 ideal run" begin
    file = "./data/GasLib-11/"
    ss = initialize_simulator(file, eos=:ideal, initial_guess_filename="")
    solver_return = run_simulator!(ss)

    @show solver_return.iterations
    @test solver_return.status == unique_physical_solution
end

@testset "test validity of compressor approximation in non ideal case" begin 
    file = "./data/GasLib-11/"
    eos = :simple_cnga

    # create and solve the potential formulation
    ss_potential = initialize_simulator(file, 
        eos = eos, 
        use_potential_formulation = true, 
        potential_ratio_coefficients = [0.0, 0.0, 0.9, 0.1], 
        initial_guess_filename = "") 
    df_potential = prepare_for_solve!(ss_potential)
    solver_return_potential = run_simulator!(ss_potential, df_potential)
    @show solver_return_potential.iterations  

    # create and solve pressure formulation 
    ss_pressure = initialize_simulator(file, 
        eos = eos,
        initial_guess_filename = "")
    df_pressure = prepare_for_solve!(ss_pressure)
    solver_return_pressure = run_simulator!(ss_pressure, df_pressure)
    @show solver_return_pressure.iterations  

    # create guess from potential solution for pressure formulation and solve
    x_guess = _create_initial_guess(ss_pressure, ss_potential)
    residual_of_guess = value!(df_pressure, x_guess)
    _check_residuals(ss_pressure, residual_of_guess)
    solver_return_pressure_with_guess = run_simulator!(ss_pressure, 
        df_pressure, x_guess = x_guess) 
    @show solver_return_pressure_with_guess.iterations 
    
    # more tests 
    err = solver_return_pressure.solution - solver_return_pressure_with_guess.solution
    @test norm(err) ≈ 0.0 atol = 1e-6 
    @test solver_return_pressure.iterations >= solver_return_pressure_with_guess.iterations
end 