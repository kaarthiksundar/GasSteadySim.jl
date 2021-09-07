@testset "test 8 node ideal run" begin
    file = "./data/8-node/"
    ss = initialize_simulator(file, eos=:ideal, initial_guess_filename="")
    solver_return = run_simulator!(ss)

    @test solver_return.status == successfull
    @test ref(ss, :node, 1, "withdrawal") * nominal_values(ss, :mass_flow) ≈ -275.00 atol = 1e-2
end

@testset "test 8 node simple CNGA run" begin
    file = "./data/8-node/"
    ss = initialize_simulator(file, eos=:simple_cnga, initial_guess_filename="")
    solver_return = run_simulator!(ss)

    @test solver_return.status == successfull
    @test ref(ss, :node, 1, "withdrawal") * nominal_values(ss, :mass_flow) ≈ -275.00 atol = 1e-2
end


@testset "test 8 node full CNGA run" begin
    file = "./data/8-node/"
    ss = initialize_simulator(file, eos=:full_cnga, initial_guess_filename="")
    solver_return = run_simulator!(ss)

    @test solver_return.status == successfull
    @test ref(ss, :node, 1, "withdrawal") * nominal_values(ss, :mass_flow) ≈ -275.00 atol = 1e-2
end
