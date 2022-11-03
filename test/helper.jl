""" helper function to check solution """
function _check_correctness(sol::Dict{String,Any}, exact_sol::Dict{String,Any}, tol = 1e-4)
    entries = ["compressor_flow", "pipe_flow", "nodal_pressure"]  
    for key in entries, i in sol[key] |> keys
        @test get(sol[key], i, 0.0) ≈ get(exact_sol[key], string(i), 0.0) rtol = 1e-2
    end
    return
end