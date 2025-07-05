


"""function assembles the residuals"""
function assemble_residual!(ss::SteadySimulator, x_dof::AbstractArray, residual_dof::AbstractArray, case::Symbol)
    _eval_junction_equations!(ss, x_dof, residual_dof)
    _eval_pipe_equations!(ss, x_dof, residual_dof, case)
    _eval_compressor_equations!(ss, x_dof, residual_dof)
    _eval_control_valve_equations!(ss, x_dof, residual_dof)
    # _eval_short_pipe_equations!(ss, x_dof, residual_dof)
    _eval_pass_through_equations!(ss, x_dof, residual_dof)
end

"""function assembles the Jacobians"""
function assemble_mat!(ss::SteadySimulator, x_dof::AbstractArray, J::AbstractArray, case::Symbol)
    fill!(J, 0)
    _eval_junction_equations_mat!(ss, x_dof, J)
    _eval_pipe_equations_mat!(ss, x_dof, J, case)
    _eval_compressor_equations_mat!(ss, x_dof, J)
    _eval_control_valve_equations_mat!(ss, x_dof, J)
    # _eval_short_pipe_equations_mat!(ss, x_dof, J)
    _eval_pass_through_equations_mat!(ss, x_dof, J)
end

# y = psi_inv(p)
function pressure_to_dof(p::Real)::Real
    return  p^3  #psi_inv(p)
end

# p = psi(y)
function dof_to_pressure(dof_val::Real)::Real
    return  cbrt(dof_val) #psi(y)
end


function  dof_to_pressure_derivative(dof_val::Real)::Real
    return  1.0/ (3 * cbrt(dof_val^2)) # psi'(y)
end

function dof_to_pressure_second_derivative(dof_val::Real)::Real
    return  -2.0 / (9 * cbrt(dof_val^5)) # psi''(y)
end

"""residual computation for junctions"""
function _eval_junction_equations!(ss::SteadySimulator, x_dof::AbstractArray, residual_dof::AbstractArray)
    @inbounds for (id, junction) in ref(ss, :node)
        eqn_no = junction["dof"]
        
        if junction["is_slack"] == 1
            pressure = junction["pressure"]
            coeff = pressure_to_dof(pressure)
            residual_dof[eqn_no] = x_dof[eqn_no] - (coeff)
        else
            r = (-junction["withdrawal"]) # inflow is positive convention
            out_edge = ref(ss, :outgoing_dofs, id)
            in_edge = ref(ss, :incoming_dofs, id)
            r -= sum(x_dof[e] for e in out_edge; init=0.0) 
            r += sum(x_dof[e] for e in in_edge; init=0.0)
            residual_dof[eqn_no] = r
        end
    end
end

# p'(x) = G(p, f) 
function G(ss::SteadySimulator, p::Real, c::Vector)::Real
    f, c0, c1, c2 = c
    rho = get_density(ss, p)
    rho_prime = get_density_prime(ss, p)
    t1 = rho^3*c2 - rho * c0 * f * abs(f)
    t2 = rho^2 - c1 * (f^2) * rho_prime
    return t1/t2
end

function dGdp(ss::SteadySimulator, p::Real, c::Vector)::Real
    f, c0, c1, c2 = c
    rho = get_density(ss, p)
    rho_prime = get_density_prime(ss, p)
    rho_double_prime = get_density_double_prime(ss, p)
    t1 = rho^3*c2 - rho * c0 * f * abs(f)
    t2 = rho^2 - c1 * (f^2) * rho_prime
    t1_prime = (3 * (rho^2) * c2 - c0 * f * abs(f)) * rho_prime
    t2_prime = 2 * rho * rho_prime - c1 * (f^2) * rho_double_prime
    return (t1_prime * t2 - t1 * t2_prime) / (t2^2)
end

function dGdf(ss::SteadySimulator, p::Real, c::Vector)::Real
    f, c0, c1, c2 = c
    rho = get_density(ss, p)
    rho_prime = get_density_prime(ss, p)
    t1 = rho^3*c2 - rho * c0 * f * abs(f)
    t2 = rho^2 - c1 * (f^2) * rho_prime
    t1_prime = - 2 * rho * c0 * abs(f)
    t2_prime = - 2  * c1 * f * rho_prime
    return (t1_prime * t2 - t1 * t2_prime) / (t2^2)  
end

# p'(x) = G(p, f) with p = psi(y) =>  y'(x) = H(y, f) =  G(psi(y), f)/ psi'(y)
# 2 pt collocation => y_fr - y_to + L/2 * ( H(y_fr, f) + H(y_to, f) ) = 0

function H(ss::SteadySimulator, y::Real, c::Vector)::Real
    p = dof_to_pressure(y)
    return G(ss, p, c)/dof_to_pressure_derivative(y)
end

function dHdy(ss::SteadySimulator, y::Real, c::Vector)::Real
    p = dof_to_pressure(y)
    return dGdp(ss, p, c) -  ( G(ss, p, c) * dof_to_pressure_second_derivative(y) ) / (dof_to_pressure_derivative(y))^2
end


function dHdf(ss::SteadySimulator, y::Real, c::Vector)::Real
    p = dof_to_pressure(y)
    return dGdf(ss, p, c)  / dof_to_pressure_derivative(y)
end

function pipe_residual(x1::Real, x2::Real)::Real
    return x1 - x2
end

function pipe_residual_derivatives(x1::Real, x2::Real)::Tuple{Real, Real}
    return 1.0, -1.0
end

"""residual computation for pipes"""
function _eval_pipe_equations!(ss::SteadySimulator, x_dof::AbstractArray, residual_dof::AbstractArray, case::Symbol)

    if case == :ode
        ode_func(y, c, x)  = H(ss, y, c)
    end
    @inbounds for (i, pipe) in ref(ss, :pipe)
        eqn_no = pipe["dof"]
        f = x_dof[eqn_no]
        fr_node = pipe["fr_node"]  
        to_node = pipe["to_node"]
        fr_dof = ref(ss, :node, fr_node, "dof")
        to_dof = ref(ss, :node, to_node, "dof")
        sin_incline = 0.065
        R1 = nominal_values(ss, :mach_num)^2 / (nominal_values(ss, :euler_num) * pipe["area"]^2)  
        # 8 degree inclination = sin(theta) approx 0.14
        # 2 degree = sin(theta) approx 0.034
        R2 = nominal_values(ss, :mach_num)^2 / (nominal_values(ss, :euler_num) *  nominal_values(ss, :froude_num)^2)
        beta = pipe["friction_factor"]  / (2 * pipe["diameter"])
        c0 = R1 * beta 
        c1 = R1 * params(ss, :inertial_bool)
        c2 = R2 * sin_incline * params(ss, :gravity_bool)
        c = [f, c0, c1, c2]
        
        if case == :ode
            prob = ODEProblem(ode_func, x_dof[fr_dof], (0, pipe["length"]), c)
            # Use ImplicitEuler() or Trapezoid() or TRBDF2()
            sol = solve(prob, Trapezoid(), save_everystep = false); 
            x_dof_to_ode = sol.u[2]
            residual_dof[eqn_no] =  pipe_residual(x_dof_to_ode, x_dof[to_dof])  # mimic p^3 term
        end

        if case ==:two_point_collocation
            residual_dof[eqn_no] =  x_dof[fr_dof]  - x_dof[to_dof]  + pipe["length"]* 0.5 * ( H(ss, x_dof[fr_dof], c) + H(ss, x_dof[to_dof], c) )
        end

    
    end
end

"""residual computation for compressor"""
function _eval_compressor_equations!(ss::SteadySimulator, x_dof::AbstractArray, residual_dof::AbstractArray)
    (!haskey(ref(ss), :compressor)) && (return)
    @inbounds for (_, comp) in ref(ss, :compressor)
        eqn_no = comp["dof"] 
        alpha = comp["c_ratio"]
        to_node = comp["to_node"]
        fr_node = comp["fr_node"]
        residual_dof[eqn_no] = pressure_to_dof(alpha) * x_dof[ref(ss, :node, fr_node, "dof")] -  x_dof[ref(ss, :node, to_node, "dof")]
    end
end

"""residual computation for control_valves"""
function _eval_control_valve_equations!(ss::SteadySimulator, x_dof::AbstractArray, residual_dof::AbstractArray)
    (!haskey(ref(ss), :control_valve)) && (return)
    @inbounds for (_, cv) in ref(ss, :control_valve)
        eqn_no = cv["dof"] 
        alpha = cv["c_ratio"]
        to_node = cv["to_node"]
        fr_node = cv["fr_node"]
        residual_dof[eqn_no] = pressure_to_dof(alpha) * x_dof[ref(ss, :node, fr_node, "dof")] -  x_dof[ref(ss, :node, to_node, "dof")]
    end
end

"""residual computation for short pipes"""
function _eval_short_pipe_equations!(ss::SteadySimulator, x_dof::AbstractArray, residual_dof::AbstractArray)
    @inbounds for (_, pipe) in get(ref(ss), :short_pipe, [])
        eqn_no = pipe["dof"]
        f = x_dof[eqn_no]
        fr_node = pipe["fr_node"]  
        to_node = pipe["to_node"]
        fr_dof = ref(ss, :node, fr_node, "dof")
        to_dof = ref(ss, :node, to_node, "dof")
        resistance = 1e-5
        residual_dof[eqn_no] = (cbrt(x_dof[fr_dof]))^2  - (cbrt(x_dof[to_dof]))^2 - f * abs(f) * resistance
    end
end

"""residual computation for pass through components"""
function _eval_pass_through_equations!(ss::SteadySimulator, x_dof::AbstractArray, residual_dof::AbstractArray)
    components = [:valve, :resistor, :loss_resistor]
    @inbounds for component in components 
        (!haskey(ref(ss), component)) && (continue)
        for (_, comp) in ref(ss, component)
            eqn_no = comp["dof"]
            fr_node = comp["fr_node"]
            to_node = comp["to_node"]
            fr_dof = ref(ss, :node, fr_node, "dof")
            to_dof = ref(ss, :node, to_node, "dof")
            residual_dof[eqn_no] = x_dof[fr_dof] - x_dof[to_dof]
        end 
    end 
end 

"""in place Jacobian computation for junctions"""
function _eval_junction_equations_mat!(ss::SteadySimulator, x_dof::AbstractArray, 
        J::AbstractArray)
    @inbounds for (id, junction) in ref(ss, :node)
        eqn_no = junction["dof"]
        
        if junction["is_slack"] == 1
            J[eqn_no, eqn_no] = 1.0
            
        else
            out_edge = ref(ss, :outgoing_dofs, id)
            in_edge = ref(ss, :incoming_dofs, id)
            for e in out_edge
                J[eqn_no, e] = -1
            end
            for e in in_edge
                J[eqn_no, e] = +1
            end
        end
    end
end


"""in place Jacobian computation for pipes"""
function _eval_pipe_equations_mat!(ss::SteadySimulator, x_dof::AbstractArray, 
        J::AbstractArray, case::Symbol)
    
    if case == :ode
        function ode_syst!(du, u, c, x)
            y, lambda_y0, lambda_f = u   
            du[1] = H(ss, y, c)  #dy
            du[2] = dHdy(ss, y, c) * lambda_y0 #dlambda_y0
            du[3] = dHdy(ss, y, c) * lambda_f + dHdf(ss, y, c) # dlambda_f
        end
    end
    @inbounds for (key, pipe) in ref(ss, :pipe)
        eqn_no = pipe["dof"] 
        f = x_dof[eqn_no]
        fr_node = pipe["fr_node"]  
        to_node = pipe["to_node"]

        eqn_fr = ref(ss, :node, fr_node, "dof")
        eqn_to = ref(ss, :node, to_node, "dof")

        sin_incline = 0.065
        R1 = nominal_values(ss, :mach_num)^2 / (nominal_values(ss, :euler_num) * pipe["area"]^2)  

        # 8 degree inclination = sin(theta) approx 0.14
        # 2 degree = sin(theta) approx 0.034
        R2 = nominal_values(ss, :mach_num)^2 / ( nominal_values(ss, :euler_num) * nominal_values(ss, :froude_num)^2 )

        beta = pipe["friction_factor"]  / (2 * pipe["diameter"])
        c0 = R1 * beta
        c1 = R1 * params(ss, :inertial_bool)
        c2 = R2 * sin_incline * params(ss, :gravity_bool)

        c = [f, c0, c1, c2]

        if case == :ode
            u0 = [x_dof[eqn_fr], 1.0, 0]
            prob = ODEProblem(ode_syst!, u0, (0, pipe["length"]), c)
            sol = solve(prob, Trapezoid(), save_everystep = false);
            y_to_ode, lambda_y0, lambda_f = sol.u[2]

            ## y
            R_x1, R_x2 = pipe_residual_derivatives(y_to_ode, x_dof[eqn_to])
            J[eqn_no, eqn_fr] = R_x1 * lambda_y0
            J[eqn_no, eqn_no] = R_x1 * lambda_f
            J[eqn_no, eqn_to] = R_x2
        end

        if case == :two_point_collocation
            J[eqn_no, eqn_fr] = 1 + pipe["length"] * 0.5 * (dHdy(ss, x_dof[eqn_fr], c))
            J[eqn_no, eqn_to] = -1 + pipe["length"] * 0.5 * (dHdy(ss, x_dof[eqn_to], c))
            J[eqn_no, eqn_no] = pipe["length"] * 0.5 * ( dHdf(ss, x_dof[eqn_to], c) + dHdf(ss, x_dof[eqn_fr], c))
        end


    end
end

"""in place Jacobian computation for compressors"""
function _eval_compressor_equations_mat!(ss::SteadySimulator, x_dof::AbstractArray, 
        J::AbstractArray)
    (!haskey(ref(ss), :compressor)) && (return)
    @inbounds for (_, comp) in ref(ss, :compressor)
        eqn_no = comp["dof"] 
        alpha = comp["c_ratio"]
        to_node = comp["to_node"]
        fr_node = comp["fr_node"]
        eqn_to = ref(ss, :node, to_node, "dof")
        eqn_fr = ref(ss, :node, fr_node, "dof")
        

        J[eqn_no, eqn_to] =   -1.0 
        J[eqn_no, eqn_fr] =  pressure_to_dof(alpha)    
    end
end

"""in place Jacobian computation for control_valves"""
function _eval_control_valve_equations_mat!(ss::SteadySimulator, x_dof::AbstractArray, 
        J::AbstractArray)
    (!haskey(ref(ss), :control_valve)) && (return)
    @inbounds for (_, cv) in ref(ss, :control_valve)
        eqn_no = cv["dof"] 
        alpha = cv["c_ratio"]
        to_node = cv["to_node"]
        fr_node = cv["fr_node"]
        eqn_to = ref(ss, :node, to_node, "dof")
        eqn_fr = ref(ss, :node, fr_node, "dof")
        J[eqn_no, eqn_to] = -1
        J[eqn_no, eqn_fr] = pressure_to_dof(alpha)  
    end
end

"""in place Jacobian computation for short pipes"""
function _eval_short_pipe_equations_mat!(ss::SteadySimulator, x_dof::AbstractArray, 
        J::AbstractArray)
    @inbounds for (_, pipe) in get(ref(ss), :short_pipe, [])
        eqn_no = pipe["dof"] 
        f = x_dof[eqn_no]
        fr_node = pipe["fr_node"]  
        to_node = pipe["to_node"]

        eqn_fr = ref(ss, :node, fr_node, "dof")
        eqn_to = ref(ss, :node, to_node, "dof")
        
        resistance = 1e-5

        J[eqn_no, eqn_fr] = 2.0 / ( 3.0 * cbrt(x_dof[eqn_fr]) )
        J[eqn_no, eqn_to] = -2.0 / ( 3.0 * cbrt(x_dof[eqn_to]) )
        J[eqn_no, eqn_no] = -2.0 * f * sign(f) * resistance
    end
end

"""in place Jacobian computation for pass through components"""
function _eval_pass_through_equations_mat!(ss::SteadySimulator, x_dof::AbstractArray, 
        J::AbstractArray)
    components = [:valve, :resistor, :loss_resistor]
    @inbounds for component in components 
        (!haskey(ref(ss), component)) && (continue)
        for (_, comp) in ref(ss, component)
            eqn_no = comp["dof"]
            to_node = comp["to_node"]
            fr_node = comp["fr_node"]
            eqn_to = ref(ss, :node, to_node, "dof")
            eqn_fr = ref(ss, :node, fr_node, "dof")

            J[eqn_no, eqn_to] = -1.0
            J[eqn_no, eqn_fr] = 1.0
        end 
    end 
end 