""" helper function to check solution """
function _check_correctness(sol::Dict{String,Any}, exact_sol::Dict{String,Any}, tol = 1e-4)
    entries = ["compressor_flow", "pipe_flow", "nodal_pressure", "control_valve_flow", "valve_flow", "resistor_flow", "loss_resistor_flow", "short_pipe_flow"]  # can add more quantities here if required
    for key in entries
        if haskey(sol, key)
            for i in keys(sol[key])
    # for key in entries, i in sol[key] |> keys
                @test sol[key][i] ≈ exact_sol[key][string(i)] rtol = tol
            end
        end
    end
    return
end

""" helper function to check all residuals (other than compressors, control_valves) are zero """ 
function _check_residuals(ss::SteadySimulator, r::Vector{Float64}, tol = 1e-4)
    components = [:node, :pipe, :valve, :resistor]
    for component in components 
        for (_, comp) in get(ref(ss), component, [])
            eqn_no = comp["dof"]
            @test r[eqn_no] ≈ 0.0 atol = tol 
        end 
    end 
end 

""" helper function to create pressure formulation initial guess from potential formulation """
function _create_initial_guess(ss_pressure::SteadySimulator, 
    ss_potential::SteadySimulator)::Vector{Float64}

    ndofs = length(ref(ss_pressure, :dof))
    guess = zeros(Float64, ndofs)

    for i in 1:length(ref(ss_potential, :dof))
        comp, id = ss_potential.ref[:dof][i]
        if comp == :node
            val = ref(ss_pressure, :is_pressure_node, id) ? ss_potential.ref[comp][id]["pressure"] : ss_potential.ref[comp][id]["potential"]
            guess[ss_pressure.ref[comp][id]["dof"]] = val
        elseif comp in [:pipe, :compressor, :valve]
            guess[ss_pressure.ref[comp][id]["dof"]] = ss_potential.ref[comp][id]["flow"]
        end 
    end 
    return guess
end 