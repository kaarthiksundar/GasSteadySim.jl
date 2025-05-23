


"""function assembles the residuals"""
function assemble_residual!(ss::SteadySimulator, x_dof::AbstractArray, residual_dof::AbstractArray)
    _eval_junction_equations!(ss, x_dof, residual_dof)
    _eval_pipe_equations!(ss, x_dof, residual_dof)
    _eval_compressor_equations!(ss, x_dof, residual_dof)
    _eval_control_valve_equations!(ss, x_dof, residual_dof)
    _eval_short_pipe_equations!(ss, x_dof, residual_dof)
    _eval_pass_through_equations!(ss, x_dof, residual_dof)
end

"""function assembles the Jacobians"""
function assemble_mat!(ss::SteadySimulator, x_dof::AbstractArray, J::AbstractArray)
    fill!(J, 0)
    _eval_junction_equations_mat!(ss, x_dof, J)
    _eval_pipe_equations_mat!(ss, x_dof, J)
    _eval_compressor_equations_mat!(ss, x_dof, J)
    _eval_control_valve_equations_mat!(ss, x_dof, J)
    _eval_short_pipe_equations_mat!(ss, x_dof, J)
    _eval_pass_through_equations_mat!(ss, x_dof, J)
end

"""residual computation for junctions"""
function _eval_junction_equations!(ss::SteadySimulator, x_dof::AbstractArray, residual_dof::AbstractArray)
    @inbounds for (id, junction) in ref(ss, :node)
        eqn_no = junction["dof"]
        
        if junction["is_slack"] == 1
            pressure = junction["pressure"]
            # coeff = ref(ss, :is_pressure_node, id) ? pressure : potential
            coeff = pressure
            residual_dof[eqn_no] = x_dof[eqn_no] - coeff
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

# p'(x) = G(p, f) => pp'(x)  = pG(p,f) => p_fr^2/2 - p_to^2/2  + (L/2) * ( p_fr*G(p_fr, f) + p_to*G(p_to, f)) = 0
# p'(x) = G(p, f) => 2pp'(x)  = 2pG(p,f)  => with y = p^2 that  y'(x) = 2 * sqrt(y)*G(sqrt(y), f)
function G(ss::SteadySimulator, p::Real, f::Real, c0::Real, c1::Real, c2::Real)::Real
    rho = get_density(ss, p)
    rho_prime = get_density_prime(ss, p)
    t1 = rho^3*c2 - rho * c0 * f * abs(f)
    t2 = rho^2 - c1 * (f^2) * rho_prime
    return t1/t2
    # return -beta * f * abs(f) / rho
end

function dGdp(ss::SteadySimulator, p::Real, f::Real, c0::Real, c1::Real, c2::Real)::Real
    rho = get_density(ss, p)
    rho_prime = get_density_prime(ss, p)
    rho_double_prime = get_density_double_prime(ss, p)
    t1 = rho^3*c2 - rho * c0 * f * abs(f)
    t2 = rho^2 - c1 * (f^2) * rho_prime
    t1_prime = (3 * (rho^2) * c2 - c0 * f * abs(f)) * rho_prime
    t2_prime = 2 * rho * rho_prime - c1 * (f^2) * rho_double_prime
    return (t1_prime * t2 - t1 * t2_prime) / (t2^2)
    # return     beta * f * abs(f) * rho_prime / (rho^2)
end

function dGdf(ss::SteadySimulator, p::Real, f::Real, c0::Real, c1::Real, c2::Real)::Real
    rho = get_density(ss, p)
    rho_prime = get_density_prime(ss, p)
    t1 = rho^3*c2 - rho * c0 * f * abs(f)
    t2 = rho^2 - c1 * (f^2) * rho_prime
    t1_prime = - 2 * rho * c0 * abs(f)
    t2_prime = - 2  * c1 * f * rho_prime
    return (t1_prime * t2 - t1 * t2_prime) / (t2^2)  
    # return -2 * beta * abs(f) / rho

end

"""residual computation for pipes"""
function _eval_pipe_equations!(ss::SteadySimulator, x_dof::AbstractArray, residual_dof::AbstractArray)
    ode_func(u, c, x)  = 2 * sqrt(abs(u)) * G(ss, sqrt(abs(u)), c[1], c[2], c[3], c[4])
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
        p_fr_sq = x_dof[fr_dof]
        p_to_sq = x_dof[to_dof]
        # prob = ODEProblem(ode_func, p_fr_sq, (0, pipe["length"]), c)
        # sol = solve(prob, Tsit5(), save_everystep = false, reltol = 1e-9, abstol = 1e-9);

        # if sol.retcode in [:MaxIters, :Unstable, :ConvergenceFailure, :Failure]
        #     println(sol.retcode)
        # end
        
        # p_to_sq_ode = sol[2]

        # residual_dof[eqn_no] =  p_to_sq_ode - p_to_sq
        # p_fr = sqrt(x_dof[fr_dof])
        # p_to = sqrt(x_dof[to_dof])
        # residual_dof[eqn_no] = x_dof[fr_dof]  - x_dof[to_dof]  - 2*pipe["length"]* beta * f * abs(f)
        # if i == 32
        #     var1 = 2*pipe["length"]* beta * f * abs(f)
        #     var2 =  pipe["length"]* 0.5 * ( 2 * p_fr * G(ss, p_fr, f, beta, c1, c2))
        #     var3 = pipe["length"]* 0.5 * ( 2 * p_to * G(ss, p_to, f, beta, c1, c2) )
        #     println("pipe $i ", var1, " ", var2, " ", var3)
        #     println("diff ",  p_fr * G(ss, p_fr, f, beta, c1, c2), " ", -beta * f * abs(f), " ", p_fr, " ", get_density(ss, p_fr))
        # end


        # residual_dof[eqn_no] =  x_dof[fr_dof]^2/2  - x_dof[to_dof]^2/2  + pipe["length"]* 0.5 * ( x_dof[fr_dof] * G(ss, x_dof[fr_dof], f, beta, c1, c2) + x_dof[to_dof] * G(ss, x_dof[to_dof], f, beta, c1, c2) )

        # p^2 as dof form
        p_fr = sqrt(x_dof[fr_dof])
        p_to = sqrt(x_dof[to_dof])
        residual_dof[eqn_no] =  x_dof[fr_dof]  - x_dof[to_dof]  + pipe["length"]* 0.5 * ( 2 * p_fr * G(ss, p_fr, f, c0, c1, c2) + 2 * p_to * G(ss, p_to, f, c0, c1, c2) )

        # pi - pj form
        # residual_dof[eqn_no] =  x_dof[fr_dof] - x_dof[to_dof]  + pipe["length"]* 0.5 * (  G(ss, x_dof[fr_dof], f, beta, c1, c2) +  G(ss, x_dof[to_dof], f, beta, c1, c2) )

        # exact form with dof p 
        # residual_dof[eqn_no] =  x_dof[fr_dof]^2/2  - x_dof[to_dof]^2/2  - pipe["length"]* beta * f * abs(f)
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
        c_0, c_1, c_2, c_3 = ss.potential_ratio_coefficients
        alpha_eff = c_0 + c_1 * alpha + c_2 * alpha^2 + c_3 * alpha^3

        is_pressure_eq = ref(ss, :is_pressure_node, fr_node) || ref(ss, :is_pressure_node, to_node)
        val = (is_pressure_eq) ? alpha : alpha_eff

        residual_dof[eqn_no] = (val^2) * x_dof[ref(ss, :node, fr_node, "dof")] -  x_dof[ref(ss, :node, to_node, "dof")]
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
        c_0, c_1, c_2, c_3 = ss.potential_ratio_coefficients
        alpha_eff = c_0 + c_1 * alpha + c_2 * alpha^2 + c_3 * alpha^3

        is_pressure_eq = ref(ss, :is_pressure_node, fr_node) || ref(ss, :is_pressure_node, to_node)
        val = (is_pressure_eq) ? alpha : alpha_eff
        residual_dof[eqn_no] = val * x_dof[ref(ss, :node, fr_node, "dof")]  - x_dof[ref(ss, :node, to_node, "dof")]
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
        residual_dof[eqn_no] = x_dof[fr_dof]^2/2  - x_dof[to_dof]^2/2 - f * abs(f) * resistance
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
            J[eqn_no, eqn_no] = 1
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
        J::AbstractArray)
    ode_func(u, c, x)  = 2 * sqrt(abs(u)) * G(ss, sqrt(abs(u)), c[1], c[2], c[3], c[4])
    @inbounds for (key, pipe) in ref(ss, :pipe)
        eqn_no = pipe["dof"] 
        f = x_dof[eqn_no]
        fr_node = pipe["fr_node"]  
        to_node = pipe["to_node"]

        eqn_fr = ref(ss, :node, fr_node, "dof")
        eqn_to = ref(ss, :node, to_node, "dof")
        p_fr = sqrt(x_dof[eqn_fr])
        p_to = sqrt(x_dof[eqn_to])

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
        p_fr_sq = x_dof[eqn_fr]
        p_to_sq = x_dof[eqn_to]
        # prob = ODEProblem(ode_func, p_fr_sq, (0, pipe["length"]), c)
        # sol = solve(prob, Tsit5(), save_everystep = false, reltol = 1e-9, abstol = 1e-9);

        # J[eqn_no, eqn_fr] = x_dof[eqn_fr] + pipe["length"] * 0.5 * (G(ss, x_dof[eqn_fr], f, beta, c1, c2) + x_dof[eqn_fr] * dGdp(ss, x_dof[eqn_fr], f, beta, c1, c2) )
        # J[eqn_no, eqn_to] = -x_dof[eqn_to] + pipe["length"] * 0.5 * (G(ss, x_dof[eqn_to], f, beta, c1, c2) + x_dof[eqn_to] * dGdp(ss, x_dof[eqn_to], f, beta, c1, c2) )
        # J[eqn_no, eqn_no] = pipe["length"] * 0.5 * (x_dof[eqn_to] * dGdf(ss, x_dof[eqn_to], f, beta, c1, c2) + x_dof[eqn_to] * dGdf(ss, x_dof[eqn_to], f, beta, c1, c2))

        # p^2 as dof
        J[eqn_no, eqn_fr] = 1 + pipe["length"] * 0.5 * (2 * G(ss, p_fr, f, c0, c1, c2) + 2 * p_fr * dGdp(ss, p_fr, f, c0, c1, c2) ) / (2 * p_fr)
        J[eqn_no, eqn_to] = -1 + pipe["length"] * 0.5 * (2 * G(ss, p_to, f, c0, c1, c2) + 2 * p_to * dGdp(ss, p_to, f, c0, c1, c2) ) / (2 * p_to)
        J[eqn_no, eqn_no] = pipe["length"] * 0.5 * ( 2 * p_to * dGdf(ss, p_to, f, c0, c1, c2) + 2 * p_fr * dGdf(ss, p_fr, f, c0, c1, c2))


        # pi - pj form with p as dof
        # J[eqn_no, eqn_fr] = 1 + pipe["length"] * 0.5 * (dGdp(ss, x_dof[eqn_fr], f, beta, c1, c2) )
        # J[eqn_no, eqn_to] = -1 + pipe["length"] * 0.5 * (dGdp(ss, x_dof[eqn_to], f, beta, c1, c2) )
        # J[eqn_no, eqn_no] = pipe["length"] * 0.5 * (dGdf(ss, x_dof[eqn_to], f, beta, c1, c2) + dGdf(ss, x_dof[eqn_to], f, beta, c1, c2))

        # exact with p^2 as dof
        # J[eqn_no, eqn_fr] =  1 #x_dof[eqn_fr]
        # J[eqn_no, eqn_to] = -1  #-x_dof[eqn_to]
        # J[eqn_no, eqn_no] = - 4  * (beta /ss.nominal_values[:euler_num]) *  pipe["length"] * abs(f)

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
        is_pressure_eq = ref(ss, :is_pressure_node, fr_node) || ref(ss, :is_pressure_node, to_node)

        c_0, c_1, c_2, c_3 = ss.potential_ratio_coefficients
        alpha_eff = c_0 + c_1 * alpha + c_2 * alpha^2 + c_3 * alpha^3

        J[eqn_no, eqn_to] = -1
        J[eqn_no, eqn_fr] = alpha^2  #is_pressure_eq ? alpha : alpha_eff
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
        is_pressure_eq = ref(ss, :is_pressure_node, fr_node) || ref(ss, :is_pressure_node, to_node)

        c_0, c_1, c_2, c_3 = ss.potential_ratio_coefficients
        alpha_eff = c_0 + c_1 * alpha + c_2 * alpha^2 + c_3 * alpha^3
        
        J[eqn_no, eqn_to] = -1
        J[eqn_no, eqn_fr] = is_pressure_eq ? alpha : alpha_eff
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

        J[eqn_no, eqn_fr] = x_dof[eqn_fr]
        J[eqn_no, eqn_to] = -x_dof[eqn_to]
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