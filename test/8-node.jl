@testset "test 8 node ideal run with LP initial guess" begin
    file = "./data/8-node/"
    ss = initialize_simulator(file, eos=:ideal)
    solver_return = run_simulator!(ss)

    @show solver_return.iterations
    @test solver_return.status == unique_physical_solution
    @test ref(ss, :node, 1, "withdrawal") * nominal_values(ss, :mass_flow) ≈ -275.00 atol = 1e-2
end

@testset "test 8 node ideal run random initial guess" begin
    file = "./data/8-node/"
    ss = initialize_simulator(file, eos=:ideal, initial_guess_from_opt=false)
    solver_return = run_simulator!(ss)

    @show solver_return.iterations
    @test solver_return.status == unique_physical_solution
    @test ref(ss, :node, 1, "withdrawal") * nominal_values(ss, :mass_flow) ≈ -275.00 atol = 1e-2
end



@testset "test 8 node simple CNGA run with LP initial guess" begin
    file = "./data/8-node/"
    ss = initialize_simulator(file, eos=:simple_cnga)
    solver_return = run_simulator!(ss)

    @show solver_return.iterations
    @test solver_return.status == unique_physical_solution
    @test ref(ss, :node, 1, "withdrawal") * nominal_values(ss, :mass_flow) ≈ -275.00 atol = 1e-2
end

@testset "test 8 node simple CNGA run with random initial guess" begin
    file = "./data/8-node/"
    ss = initialize_simulator(file, eos=:simple_cnga, initial_guess_from_opt=false)
    solver_return = run_simulator!(ss)

    @show solver_return.iterations
    @test solver_return.status == unique_physical_solution
    @test ref(ss, :node, 1, "withdrawal") * nominal_values(ss, :mass_flow) ≈ -275.00 atol = 1e-2
end


@testset "test 8 node full CNGA run with LP initial guess" begin
    file = "./data/8-node/"
    ss = initialize_simulator(file, eos=:full_cnga)
    solver_return = run_simulator!(ss)

    @show solver_return.iterations
    @test solver_return.status == unique_physical_solution
    @test ref(ss, :node, 1, "withdrawal") * nominal_values(ss, :mass_flow) ≈ -275.00 atol = 1e-2
end

@testset "test 8 node full CNGA run with random initial guess" begin
    file = "./data/8-node/"
    ss = initialize_simulator(file, eos=:full_cnga, initial_guess_from_opt=false)
    solver_return = run_simulator!(ss)

    @show solver_return.iterations
    @test solver_return.status == unique_physical_solution
    @test ref(ss, :node, 1, "withdrawal") * nominal_values(ss, :mass_flow) ≈ -275.00 atol = 1e-2
end
