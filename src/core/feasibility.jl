function add_flow_bounds_to_ref!(ss::SteadySimulator)
    for (i, pipe) in ref(ss, :pipe)
        fr_node = pipe["fr_node"]
        to_node = pipe["to_node"]
        fr_node_p_min = ref(ss, :node, fr_node, "min_pressure")
        fr_node_p_max = ref(ss, :node, fr_node, "max_pressure")
        to_node_p_min = ref(ss, :node, to_node, "min_pressure")
        to_node_p_max = ref(ss, :node, to_node, "max_pressure")
        c = nominal_values(ss, :mach_num)^2 / nominal_values(ss, :euler_num) 
        b1, b2 = get_eos_coeffs(ss)
        resistance = pipe["friction_factor"] * pipe["length"] * c / (2 * pipe["diameter"] * pipe["area"]^2)
    
        beta = 1/resistance
        if isinf(beta) 
            beta = 1e5
        end
        p_sqr_max = fr_node_p_max^2 - to_node_p_min^2 
        p_cube_max = fr_node_p_max^3 - to_node_p_min^3 

        p_sqr_min = to_node_p_max^2 - fr_node_p_min^2 
        p_cube_min = to_node_p_max^3 - fr_node_p_min^3 

        pipe["max_flow"] = sqrt(beta * ((b1/2) * p_sqr_max + (b2/3) * p_cube_max))
        pipe["min_flow"] = -sqrt(beta * ((b1/2) * p_sqr_min + (b2/3) * p_cube_min))
    end 
end 

function construct_feasibility_model!(ss::SteadySimulator; feasibility_model::Symbol=:milp, num_partitions::Int=4)
    milp = (feasibility_model == :milp) ? true : false
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
    var[:w] = @variable(m, [i in keys(ref(ss, :node))], 
        base_name = "withdrawal")
    var[:f_pipe] = @variable(m, [i in keys(ref(ss, :pipe))],
        lower_bound = ref(ss, :pipe, i, "min_flow"),
        upper_bound = ref(ss, :pipe, i, "max_flow"),
        base_name = "f_pipe") 
    var[:f_compressor] = @variable(m, [i in keys(ref(ss, :compressor))], 
        lower_bound = 0.0,
        base_name = "f_compressor") 
    var[:f_control_valve] = @variable(m, [i in keys(ref(ss, :control_valve))], 
        lower_bound = 0.0,
        base_name = "f_compressor") 

    # auxiliary variables 
    if (b2 == 0)
        var[:pi] = @variable(m, [i in keys(ref(ss, :node))], 
            lower_bound = (b1/2) * ref(ss, :node, i, "min_pressure")^2, 
            base_name = "pi_ideal")
    else 
        var[:pi] = @variable(m, [i in keys(ref(ss, :node))], 
            lower_bound = (b1/2) * ref(ss, :node, i, "min_pressure")^2 + (b2/3) * ref(ss, :node, i, "min_pressure")^3, 
            base_name = "p_non_ideal")
    end 
    var[:f_abs_f] = @variable(m, [i in keys(ref(ss, :pipe))], base_name = "f_abs_f")

    # relaxation constraints 
    for (i, node) in ref(ss, :node)
        _, val = control(ss, :node, i)
        if node["is_slack"] == 1
            rhs = (b1/2) * val^2 + (b2/3) * val^3
            @constraint(m, var[:pi][i] == rhs) 
            continue
        end 
        min_pressure = ref(ss, :node, i, "min_pressure") 
        max_pressure = ref(ss, :node, i, "max_pressure")
        partition = collect(range(min_pressure, max_pressure, length = num_partitions))
        if (b2 == 0)
            f = p -> (b1/2) * p^2
            f_dash = p -> b1 * p
            construct_univariate_relaxation!(m, f, var[:p][i], var[:pi][i], partition, milp; f_dash=f_dash)
        else 
            f = p -> (b1/2) * p^2 + (b2/3) * p^3
            f_dash = p -> b1 * p + b2 * p^2
            construct_univariate_relaxation!(m, f, var[:p][i], var[:pi][i], partition, milp; f_dash=f_dash)
        end
    end 

    for (i, pipe) in ref(ss, :pipe)
        min_flow = pipe["min_flow"]
        max_flow = pipe["max_flow"]
        partition = []
        if (min_flow < 0)
            partition = collect(range(min_flow, 0.0, length = Int(num_partitions/2)))
            append!(partition, collect(range(0.0, max_flow, length = Int(num_partitions/2)))[2:end])
        else 
            partition = collect(range(min_flow, max_flow, length = num_partitions))
        end 
        construct_univariate_relaxation!(m, a->a*abs(a), var[:f_pipe][i], var[:f_abs_f][i], partition, milp; 
            f_dash=a->2*a*sign(a))
    end 

    # nodal constraints 
    con[:node] = Dict{Int,Any}()
    con[:node_inputs] = Dict{Int,Any}()
    for (i, node) in ref(ss, :node)
        _, val = control(ss, :node, i)
        if node["is_slack"] == 1 
            con[:node_inputs][i] = @constraint(m, var[:p][i] == val) 
        else 
            con[:node_inputs][i] = @constraint(m, var[:w][i] == val)
        end 
        inflow = 0.0 
        outflow = var[:w][i]
        for k in ref(ss, :incoming_pipes, i)
            inflow += var[:f_pipe][k]
        end 
        for k in ref(ss, :incoming_compressors, i)
            inflow += var[:f_compressor][k]
        end 
        for k in ref(ss, :incoming_control_valves, i)
            inflow += var[:f_control_valve][k]
        end 
        for k in ref(ss, :outgoing_pipes, i)
            outflow += var[:f_pipe][k]
        end 
        for k in ref(ss, :outgoing_compressors, i)
            outflow += var[:f_compressor][k]
        end 
        for k in ref(ss, :outgoing_control_valves, i)
            outflow += var[:f_control_valve][k]
        end 
        con[:node][i] = @constraint(m, inflow - outflow == 0)
    end 

    # balance constraint
    con[:balance] = @constraint(m, sum(var[:w]) == 0)

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

    # control valve constraints 
    if haskey(ref(ss), :control_valve)
        con[:control_valve] = Dict{Int,Any}()
        for (i, cv) in ref(ss, :control_valve)
            ctrl, val = control(ss, :control_valve, i)
            if ctrl  == c_ratio_control
                to_node = cv["to_node"]
                fr_node = cv["fr_node"]
                con[:control_valve][i] = @constraint(m, var[:p][to_node] == val * var[:p][fr_node])
            elseif ctr == flow_control
                con[:control_valve][i] = @constraint(m, var[:f_control_valve][i] == val)
            elseif ctr == discharge_pressure_control
                to_node = cv["to_node"]
                con[:control_valve][i] = @constraint(m, var[:p][to_node] == val)
            end
        end 
    end

    # pipe constraints 
    con[:pipe] = Dict{Int,Any}()
    for (i, pipe) in ref(ss, :pipe)
        fr_node = pipe["fr_node"]  
        to_node = pipe["to_node"]
        c = nominal_values(ss, :mach_num)^2 / nominal_values(ss, :euler_num) 
        resistance = pipe["friction_factor"] * pipe["length"] * c / (2 * pipe["diameter"] * pipe["area"]^2)
        con[:pipe][i] = @constraint(m, var[:pi][fr_node] - var[:pi][to_node] - var[:f_abs_f][i] * resistance == 0)
    end 

end