@testset "test GasLib-11 ideal run" begin
    file = "./data/GasLib-11/"
    ss = initialize_simulator(file, eos=:ideal, initial_guess_filename="")
    solver_return = run_simulator!(ss)

    @show solver_return.stats
    @test solver_return.status == physical_solution
end

@testset "test GasLib-11 non ideal case" begin 
    file = "./data/GasLib-11/"
    ss = initialize_simulator(file, eos=:simple_cnga, initial_guess_filename="")
    solver_return = run_simulator!(ss)

    @show solver_return.stats
    @test solver_return.status == physical_solution 
end