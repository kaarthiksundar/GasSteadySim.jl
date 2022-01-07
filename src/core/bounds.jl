function _add_flow_bounds_to_ref!(ss::SteadySimulator)
    for (_, pipe) in ref(ss, :pipe)
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