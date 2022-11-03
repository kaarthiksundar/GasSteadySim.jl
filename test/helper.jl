""" helper function to check solution """
function _check_correctness(sol::Dict{String,Any}, exact_sol::Dict{String,Any}, tol = 1e-4)
    entries = ["compressor_flow", "pipe_flow", "nodal_pressure"]  
    for key in entries, i in sol[key] |> keys
        @test sol[key][i] ≈ exact_sol[key][string(i)] rtol = tol
    end
    return
end