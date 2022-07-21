function set_bounds(x::JuMP.VariableRef, lb::Number, ub::Number)
    JuMP.set_lower_bound(x, lb)
    JuMP.set_upper_bound(x, ub)
end 

function populate_lp_model!(ss::SteadySimulator; 
    num_pressure_partitions::Int64 = 5, 
    num_positive_flow_partitions::Int64 = 5)
    model = ss.lp_relax
    m = model.model
    var = model.variables 
    sol = model.sol 
    residual = model.residuals

    # initialize variables (x contains [potentials, lifted flows for pipes, flows for the rest of the components])
    var[:x] = x = JuMP.@variable(m, [i in keys(ref(ss, :dof))])
    non_slack_node_dofs = filter(
        tuple -> first(last(tuple)) == :node && !ref(ss, :node, last(last(tuple)), "is_slack"), 
        ref(ss, :dof)
    )
    pressure_node_dofs = filter(
        tuple -> first(last(tuple)) == :node && ref(ss, :is_pressure_node, last(last(tuple))), 
        ref(ss, :dof)
    )
    pipe_dofs = filter(tuple -> first(last(tuple)) == :pipe, ref(ss, :dof))
    var[:pressure] = pressure = JuMP.@variable(m, [i in keys(pressure_node_dofs)])
    val[:flow] = flow = JuMP.@variable(m, [i in keys(pipe_dofs)])
    var[:t] = t = JuMP.@variable(m, [i in keys(non_slack_node_dofs)], lower_bound = 0.0)

    # set bounds
    for (dof, val) in ref(ss, :dof)
        id = val[2]
        if val[1] == :node 
            lb = get_potential(ss, ref(ss, :node, id, "min_pressure"))
            ub = get_potential(ss, ref(ss, :node, id, "max_pressure"))
            set_bounds(x[dof], lb, ub)
            if ref(ss, :is_pressure_node, id)
                set_bounds(pressure[dof], ref(ss, :node, id, "min_pressure"), ref(ss, :node, id, "max_pressure"))
            end 
        elseif val[1] == :pipe 
            set_bounds(flow[dof], ref(ss, :pipe, id, "min_flow"), ref(ss, :pipe, id, "max_flow"))
        else 
            continue 
        end 
    end 

    # add node constraints 
    for (node_id, node) in ref(ss, :node)
        eqn_no = node["dof"]
        ctrl_type, val = control(ss, :node, node_id) # val is withdrawal or pressure

        if ctrl_type == pressure_control
            potential = get_potential(ss, val)
            JuMP.@constraint(m, x[eqn_no] == potential)
            if ref(ss, :is_pressure_node, node_id)
                JuMP.@constraint(m, pressure[eqn_no] == val)
            end 
        end

        if  ctrl_type == flow_control
            out_edge = ref(ss, :outgoing_dofs, node_id)
            in_edge = ref(ss, :incoming_dofs, node_id)
            balance_expr = sum(flow[e] for e in in_edge; init=0.0) - sum(flow[e] for e in out_edge; init=0.0) - val
            JuMP.@constraint(m, balance_expr >= -t[eqn_no])
            JuMP.@constraint(m, balance_expr <= t[eqn_no])
        end

        if ref(ss, :is_pressure_node, node_id)
            lb = JuMP.get_lower_bound(pressure[eqn_no])
            ub = JuMP.get_upper_bound(pressure[eqn_no])
            partition = collect(range(start=lb, stop=ub, length=num_pressure_partitions+1))
            construct_univariate_relaxation!(m, 
                p -> get_potential(ss, p), 
                pressure[eqn_no], x[eqn_no], 
                partition, false, 
                f_dash=p->get_potential_derivative(ss, p)
            )
        end 
    end


    # add pipe constraints
    for (_, pipe) in ref(ss, :pipe)
        eqn_no = pipe["dof"]
        f = x[eqn_no]
        fr_node = pipe["fr_node"]  
        to_node = pipe["to_node"]
        fr_dof = ref(ss, :node, fr_node, "dof")
        to_dof = ref(ss, :node, to_node, "dof")
        is_fr_pressure_node = ref(ss, :is_pressure_node, fr_node)
        is_to_pressure_node = ref(ss, :is_pressure_node, to_node)
        pi_fr = (is_fr_pressure_node) ? get_potential(ss, x_dof[fr_dof]) : x_dof[fr_dof] 
        pi_to = (is_to_pressure_node) ? get_potential(ss, x_dof[to_dof]) : x_dof[to_dof] 
        c = nominal_values(ss, :mach_num)^2 / nominal_values(ss, :euler_num) 

        resistance = pipe["friction_factor"] * pipe["length"] * c / (2 * pipe["diameter"] * pipe["area"]^2)
        residual_dof[eqn_no] = pi_fr - pi_to - f * abs(f) * resistance
    end
    

end 

