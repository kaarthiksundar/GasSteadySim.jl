

function _check_correctness(sol::Dict{String,Any}, exact_sol::Dict{String,Any})

    entries = ["compressor_flow", "pipe_flow", "nodal_pressure"]  # can add more quantities here if required
    for key in entries
        for i in keys(sol[key])
            # @show sol[key][i] - exact_sol[key][string(i)]
            @test isapprox(sol[key][i], exact_sol[key][string(i)]; rtol = 5e-1)
        end
    end
    return
end