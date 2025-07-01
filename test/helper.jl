""" helper function to check solution """
function _check_correctness(sol::Dict{String,Any}, exact_sol::Dict{String,Any}, tol = 1e-3)
    entries = ["compressor_flow", "pipe_flow", "nodal_pressure"]  
    for key in entries, i in sol[key] |> keys
        @test sol[key][i] ≈ exact_sol[key][string(i)] rtol = tol
    end
    return
end

""" helper function to check all residuals (other than compressors, control_valves) are zero """ 
function _check_residuals(ss::SteadySimulator, r::Vector{Float64}, tol = 1e-3)
    components = [:node, :pipe, :valve, :resistor]
    for component in components 
        for (_, comp) in get(ref(ss), component, [])
            eqn_no = comp["dof"]
            @test r[eqn_no] ≈ 0.0 atol = tol 
        end 
    end 
end 

