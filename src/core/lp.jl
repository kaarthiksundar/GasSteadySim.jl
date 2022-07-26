function _set_bounds(x::JuMP.VariableRef, lb::Number, ub::Number)
    JuMP.set_lower_bound(x, lb)
    JuMP.set_upper_bound(x, ub)
end 

function _populate_lp_model!(ss::SteadySimulator; 
    num_pressure_partitions::Int64 = 5, 
    num_positive_flow_partitions::Int64 = 5)
    model = ss.lp_relax
    m = model.model
    var = model.variables 

    # initialize variables (x contains [potentials, lifted flows for pipes, flows for the rest of the components])
    var[:x] = x = JuMP.@variable(m, [i in keys(ref(ss, :dof))])
    non_slack_node_dofs = filter(
        tuple -> first(last(tuple)) == :node && ref(ss, :node, last(last(tuple)), "is_slack") == false, 
        ref(ss, :dof)
    )
    pressure_node_dofs = filter(
        tuple -> first(last(tuple)) == :node && ref(ss, :is_pressure_node, last(last(tuple))), 
        ref(ss, :dof)
    )
    pipe_dofs = filter(tuple -> first(last(tuple)) == :pipe, ref(ss, :dof))
    var[:pressure] = pressure = JuMP.@variable(m, [i in keys(pressure_node_dofs)])
    var[:flow] = flow = JuMP.@variable(m, [i in keys(pipe_dofs)])
    var[:t] = t = JuMP.@variable(m, [i in keys(non_slack_node_dofs)], lower_bound = 0.0)

    # set bounds
    for (dof, val) in ref(ss, :dof)
        id = val[2]
        if val[1] == :node 
            lb = get_potential(ss, ref(ss, :node, id, "min_pressure"))
            ub = get_potential(ss, ref(ss, :node, id, "max_pressure"))
            _set_bounds(x[dof], lb, ub)
            if ref(ss, :is_pressure_node, id)
                _set_bounds(pressure[dof], ref(ss, :node, id, "min_pressure"), ref(ss, :node, id, "max_pressure"))
            end 
        elseif val[1] == :pipe 
            _set_bounds(flow[dof], ref(ss, :pipe, id, "min_flow"), ref(ss, :pipe, id, "max_flow"))
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
            in_expr = sum((first(ref(ss, :dof, e)) == :pipe) ? flow[e] : x[e] for e in in_edge; init=0.0)
            out_expr = sum((first(ref(ss, :dof, e)) == :pipe) ? flow[e] : x[e] for e in out_edge; init=0.0)
            balance_expr = in_expr - out_expr - val
            JuMP.@constraint(m, balance_expr >= -t[eqn_no])
            JuMP.@constraint(m, balance_expr <= t[eqn_no])
        end

        if ref(ss, :is_pressure_node, node_id)
            lb = JuMP.lower_bound(pressure[eqn_no])
            ub = JuMP.upper_bound(pressure[eqn_no])
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
        pi_fr = x[fr_dof]
        pi_to = x[to_dof] 
        c = nominal_values(ss, :mach_num)^2 / nominal_values(ss, :euler_num) 
        resistance = pipe["friction_factor"] * pipe["length"] * c / (2 * pipe["diameter"] * pipe["area"]^2)
        JuMP.@constraint(m, pi_fr - pi_to == f * resistance)
        
        lb = JuMP.lower_bound(flow[eqn_no])
        ub = JuMP.upper_bound(flow[eqn_no])
        partition = [
            collect(range(start=lb, stop=0.0, length=num_positive_flow_partitions+1))..., 
            collect(range(start=0.0, stop=ub, length=num_positive_flow_partitions+1))[2:end]...
        ]
        construct_univariate_relaxation!(m, 
            y -> y * abs(y), 
            flow[eqn_no], f, 
            partition, false, 
            f_dash=y->2 * abs(y)
        )
    end

    # add compressor constraints 
    for (comp_id, comp) in get(ref(ss), :compressor, [])
        _ = comp["dof"] 
        _, cmpr_val = control(ss, :compressor, comp_id)
        to_node = comp["to_node"]
        fr_node = comp["fr_node"]
        
        is_pressure_eq = ref(ss, :is_pressure_node, fr_node) || ref(ss, :is_pressure_node, to_node)
        val = (is_pressure_eq) ? cmpr_val : cmpr_val^2
        if (is_pressure_eq)
            JuMP.@constraint(m, pressure[ref(ss, :node, fr_node, "dof")] * val == pressure[ref(ss, :node, to_node, "dof")])
        else 
            JuMP.@constraint(m, x[ref(ss, :node, fr_node, "dof")] * val == x[ref(ss, :node, to_node, "dof")])
        end 
    end

    # add control valve constraints 
    for (cv_id, cv) in get(ref(ss), :control_valve, [])
        _ = cv["dof"] 
        _, cv_val = control(ss, :control_valve, cv_id)
        to_node = cv["to_node"]
        fr_node = cv["fr_node"]
        
        is_pressure_eq = ref(ss, :is_pressure_node, fr_node) || ref(ss, :is_pressure_node, to_node)
        val = (is_pressure_eq) ? cv_val : cv_val^2
        if (is_pressure_eq)
            JuMP.@constraint(m, val * pressure[ref(ss, :node, fr_node, "dof")] == pressure[ref(ss, :node, to_node, "dof")])
        else 
            JuMP.@constraint(m, x[ref(ss, :node, fr_node, "dof")] * val == x[ref(ss, :node, to_node, "dof")])
        end 
    end

    # add pass-through component constraints
    components = [:valve, :resistor, :loss_resistor, :short_pipe]
    for component in components 
        (!haskey(ref(ss), component)) && (continue)
        for (_, comp) in ref(ss, component)
            _ = comp["dof"]
            fr_node = comp["fr_node"]
            to_node = comp["to_node"]
            fr_dof = ref(ss, :node, fr_node, "dof")
            to_dof = ref(ss, :node, to_node, "dof")
            is_fr_pressure_node = ref(ss, :is_pressure_node, fr_node)
            is_to_pressure_node = ref(ss, :is_pressure_node, to_node)
            if (is_fr_pressure_node && is_to_pressure_node)
                JuMP.@constraint(m, pressure[fr_dof] == pressure[to_dof])
            else
                JuMP.@constraint(m, x[fr_dof] == x[to_dof])
            end 
        end 
    end 

    # add objective 
    JuMP.@objective(m, Min, sum(t))
end 

# solve the LP relaxation model
function _solve_lp_model!(ss::SteadySimulator)
    JuMP.optimize!(ss.lp_relax.model)
end

_is_optimal(ss::SteadySimulator) = JuMP.termination_status(ss.lp_relax.model) == MOI.OPTIMAL

# populate the LP solution for an initial guess with the residuals
function _populate_lp_solution!(ss::SteadySimulator)
    om = ss.lp_relax
    m = om.model
    (JuMP.termination_status(m) != MOI.OPTIMAL) && (return)
    sol(om)[:x] = Dict{Int, Any}()

    for (node_id, node) in ref(ss, :node)
        eqn_no = node["dof"]
        sol(om, :x)[eqn_no] = (ref(ss, :is_pressure_node, node_id)) ? JuMP.value(var(om, :pressure, eqn_no)) : JuMP.value(var(om, :x, eqn_no))
    end
    
    for (_, pipe) in ref(ss, :pipe)
        eqn_no = pipe["dof"]
        sol(om, :x)[eqn_no] = JuMP.value(var(om, :flow, eqn_no))
    end 
    
    components = [:compressor, :control_valve, :valve, :resistor, :loss_resistor, :short_pipe]
    for component in components 
        (!haskey(ref(ss), component)) && (continue)
        for (_, comp) in ref(ss, component)
            eqn_no = comp["dof"]
            sol(om, :x)[eqn_no] = JuMP.value(var(om, :x, eqn_no))
        end 
    end
    
end

