function add_flow_bounds_to_ref!(ss::SteadySimulator)
    for (_, pipe) in ref(ss, :pipe)
        fr_node = pipe["fr_node"]
        to_node = pipe["to_node"]
        fr_node_p_min = ref(ss, :node, fr_node, "min_pressure")
        fr_node_p_max = ref(ss, :node, fr_node, "max_pressure")
        to_node_p_min = ref(ss, :node, to_node, "min_pressure")
        to_node_p_max = ref(ss, :node, to_node, "max_pressure")
        c = ss.nominal_values[:mass_flow]^2  / (ss.nominal_values[:pressure] * ss.nominal_values[:density])

        b1, b2 = get_eos_coeffs(ss)
        resistance = pipe["friction_factor"] * pipe["length"] * c / (2 * pipe["diameter"] * pipe["area"]^2)
        beta = 1/resistance
        p_sqr_max = fr_node_p_max^2 - to_node_p_min^2 
        p_cube_max = fr_node_p_max^3 - to_node_p_min^3 

        p_sqr_min = to_node_p_max^2 - fr_node_p_min^2 
        p_cube_min = to_node_p_max^3 - fr_node_p_min^3 

        pipe["max_flow"] = sqrt(beta * ((b1/2) * p_sqr_max + (b2/3) * p_cube_max))
        pipe["min_flow"] = -sqrt(beta * ((b1/2) * p_sqr_min + (b2/3) * p_cube_min))
    end 
end 

function construct_feasibility_model!(ss::SteadySimulator)
    add_flow_bounds_to_ref!(ss)
    var = ss.variables
    con = ss.constraints
    m = ss.feasibility_model
    b1, b2 = get_eos_coeffs(ss)

    # state variables
    var[:p] = @variable(m, [i in keys(ref(ss, :node))], 
        lower_bound = ref(ss, :node, i, "min_pressure"), 
        upper_bound = ref(ss, :node, i, "max_pressure"), 
        base_name = "p") 
    var[:f_pipe] = @variable(m, [i in keys(ref(ss, :pipe))],
        lower_bound = ref(ss, :pipe, i, "min_flow"),
        upper_bound = ref(ss, :pipe, i, "max_flow"),
        base_name = "f_pipe") 
    var[:f_compressor] = @variable(m, [i in keys(ref(ss, :compressor))], 
        lower_bound = 0.0,
        base_name = "f_compressor") 

    # auxiliary variables 
    var[:p_sqr] = @variable(m, [i in keys(ref(ss, :node))], base_name = "p_sqr")
    (b2 != 0.0) && (var[:p_cube] = @variable(m, [i in keys(ref(ss, :node))], base_name = "p_cube"))
    var[:f_abs_f] = @variable(m, [i in keys(ref(ss, :pipe))], base_name = "f_abs_f")

    # relaxation constraints 
    for (i, node) in ref(ss, :node)
        _, val = control(ss, :node, i)
        if node["is_slack"] == 1
            @constraint(m, var[:p_sqr][i] == val^2) 
            if (b2 != 0)
                @constraint(m, var[:p_cube][i] == val^3) 
            end 
            continue
        end 
        min_pressure = ref(ss, :node, i, "min_pressure") 
        max_pressure = ref(ss, :node, i, "max_pressure")
        partition = collect(range(min_pressure, max_pressure, length = 6))
        construct_univariate_relaxation!(m, a->a^2, var[:p][i], var[:p_sqr][i], partition, false)
        if (b2 != 0)
            construct_univariate_relaxation!(m, a->a^3, var[:p][i], var[:p_cube][i], partition, false)
        end
    end 

    for (i, pipe) in ref(ss, :pipe)
        min_flow = pipe["min_flow"]
        max_flow = pipe["max_flow"]
        partition = []
        if (min_flow < 0)
            partition = collect(range(min_flow, 0.0, length = 3))
            append!(partition, collect(range(0.0, max_flow, length = 3))[2:end])
        else 
            partition = collect(range(min_flow, max_flow, length = 6))
        end 
        construct_univariate_relaxation!(m, a->a*abs(a), var[:f_pipe][i], var[:f_abs_f][i], partition, false; 
            f_dash=a->2*a*sign(a))
    end 


    # nodal constraints 
    con[:node] = Dict{Int,Any}()
    for (i, node) in ref(ss, :node)
        _, val = control(ss, :node, i)
        if node["is_slack"] == 1 
            con[:node][i] = @constraint(m, var[:p][i] == val) 
            continue 
        end
        inflow = 0.0 
        outflow = val
        for k in ref(ss, :incoming_pipes, i)
            inflow += var[:f_pipe][k]
        end 
        for k in ref(ss, :incoming_compressors, i)
            inflow += var[:f_compressor][k]
        end 
        for k in ref(ss, :outgoing_pipes, i)
            outflow += var[:f_pipe][k]
        end 
        for k in ref(ss, :outgoing_compressors, i)
            outflow += var[:f_compressor][k]
        end 
        con[:node][i] = @constraint(m, inflow - outflow == 0)
    end 

    # compressor constraints 
    con[:compressor] = Dict{Int,Any}()
    for (i, comp) in ref(ss, :compressor)
        ctrl, val = control(ss, :compressor, i)
        if ctrl  == c_ratio_control
            to_node = comp["to_node"]
            fr_node = comp["fr_node"]
            con[:compressor][i] = @constraint(m, var[:p][to_node] == val * var[:p][fr_node])
        elseif ctr == flow_control
            con[:compressor][i] = @constraint(m, var[:f_compressor][i] == val)
        elseif ctr == discharge_pressure_control
            to_node = comp["to_node"]
            con[:compressor][i] = @constraint(m, var[:p][to_node] == val)
        end
    end 

    # pipe constraints 
    con[:pipe] = Dict{Int,Any}()
    for (i, pipe) in ref(ss, :pipe)
        fr_node = pipe["fr_node"]  
        to_node = pipe["to_node"]
        b1, b2 = get_eos_coeffs(ss)
        c = ss.nominal_values[:mass_flow]^2  / (ss.nominal_values[:pressure] * ss.nominal_values[:density])
        resistance = pipe["friction_factor"] * pipe["length"] * c / (2 * pipe["diameter"] * pipe["area"]^2)
        sqr_term = (b1/2) * (var[:p_sqr][fr_node] - var[:p_sqr][to_node])

        if (b2 != 0)
            cube_term = (b2/3) * (var[:p_cube][to_node] - var[:p_cube][to_node])
            con[:pipe][i] = @constraint(m,  sqr_term + cube_term - var[:f_abs_f][i] * resistance == 0)
        else 
            con[:pipe][i] = @constraint(m,  sqr_term - var[:f_abs_f][i] * resistance == 0)
        end 
    end 

end