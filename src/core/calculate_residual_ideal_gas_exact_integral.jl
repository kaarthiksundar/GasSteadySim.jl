

function get_applicable_residual(ss::SteadySimulator, x_dof::AbstractArray)::Tuple{Real, Real}

    if ss.params[:inertial_bool] == false && ss.params[:gravity_bool] == false
        return pipe_equations_no_gravity_no_inertia(ss, x_dof)
    end

    
    if ss.params[:inertial_bool] == true && ss.params[:gravity_bool] == false
        return pipe_equations_no_gravity_with_inertia(ss, x_dof)
    end


    if ss.params[:inertial_bool] == false && ss.params[:gravity_bool] == true
        return pipe_equations_with_gravity_no_inertia(ss, x_dof)
    end

    
    if ss.params[:inertial_bool] == true && ss.params[:gravity_bool] == true
        return pipe_equations_with_gravity_with_inertia(ss, x_dof)
    end
end


"""residual computation for pipes"""
function pipe_equations_no_gravity_no_inertia(ss::SteadySimulator, x_dof::AbstractArray)::Tuple{Real, Real}
    err = Vector{Float64}()
    @inbounds for (_, pipe) in ref(ss, :pipe)
        eqn_no = pipe["dof"]
        f = x_dof[eqn_no]
        fr_node = pipe["fr_node"]  
        to_node = pipe["to_node"]
        fr_dof = ref(ss, :node, fr_node, "dof")
        to_dof = ref(ss, :node, to_node, "dof")

        R1 = nominal_values(ss, :mach_num)^2 / (nominal_values(ss, :euler_num) * pipe["area"]^2)  
        beta = pipe["friction_factor"]  / (2 * pipe["diameter"])
        R1_bar = R1 / nominal_values(ss, :euler_num)
        
        p_fr = dof_to_pressure(x_dof[fr_dof])
        p_to = dof_to_pressure(x_dof[to_dof])
        
        var = abs(p_fr^2/2  - p_to^2/2  -  pipe["length"] * beta * (R1_bar) * f * abs(f) )
        push!(err, var)
    end
    err_max  = maximum(err)
    err_rms = sqrt(sum(err.^2)/ length(err))
    return err_max, err_rms
end

function pipe_equations_no_gravity_with_inertia(ss::SteadySimulator, x_dof::AbstractArray)::Tuple{Real, Real}
    err = Vector{Float64}()
    @inbounds for (_, pipe) in ref(ss, :pipe)
        eqn_no = pipe["dof"]
        f = x_dof[eqn_no]
        fr_node = pipe["fr_node"]  
        to_node = pipe["to_node"]
        fr_dof = ref(ss, :node, fr_node, "dof")
        to_dof = ref(ss, :node, to_node, "dof")

        
        R1 = nominal_values(ss, :mach_num)^2 / (nominal_values(ss, :euler_num) * pipe["area"]^2)  
 
        beta = pipe["friction_factor"]  / (2 * pipe["diameter"])
        R1_bar = R1 / nominal_values(ss, :euler_num)

        p_fr = dof_to_pressure(x_dof[fr_dof])
        p_to = dof_to_pressure(x_dof[to_dof])

        var =  abs(p_fr^2/2  - p_to^2/2   - (R1_bar) * (f^2)  * log( abs(p_fr / p_to) )- pipe["length"] * beta * (R1_bar) * f * abs(f))
        push!(err, var)
    end
    err_max  = maximum(err)
    err_rms = sqrt(sum(err.^2)/ length(err))
    return err_max, err_rms
end

function pipe_equations_with_gravity_no_inertia(ss::SteadySimulator, x_dof::AbstractArray)::Tuple{Real, Real}
    err = Vector{Float64}()
    @inbounds for (_, pipe) in ref(ss, :pipe)
        eqn_no = pipe["dof"]
        f = x_dof[eqn_no]
        fr_node = pipe["fr_node"]  
        to_node = pipe["to_node"]
        fr_dof = ref(ss, :node, fr_node, "dof")
        to_dof = ref(ss, :node, to_node, "dof")

        R1 = nominal_values(ss, :mach_num)^2 / (nominal_values(ss, :euler_num) * pipe["area"]^2)  
        R2 = nominal_values(ss, :mach_num)^2 / ( nominal_values(ss, :euler_num) * nominal_values(ss, :froude_num)^2 )
        beta = pipe["friction_factor"]  / (2 * pipe["diameter"])

        R1_bar = R1 / nominal_values(ss, :euler_num)
        R2_bar = R2 * nominal_values(ss, :euler_num)
        # 8 degree inclination = sin(theta) approx 0.14
        # 2 degree = sin(theta) approx 0.034
        sin_incline = 0.065
        p_fr = dof_to_pressure(x_dof[fr_dof])
        p_to = dof_to_pressure(x_dof[to_dof])
        gamma = 2 * sin_incline * R2_bar * pipe["length"]
        if gamma < 1e-9
            gravity_factor = 1
        else
            gravity_factor = (exp(gamma) - 1) /(gamma)
        end

        var =  abs( exp(gamma) * p_fr^2/2  - p_to^2/2  - pipe["length"] * beta * R1_bar * f * abs(f) * gravity_factor )
        push!(err, var)
    end
    err_max  = maximum(err)
    err_rms = sqrt(sum(err.^2)/ length(err))
    return err_max, err_rms
end

function pipe_equations_with_gravity_with_inertia(ss::SteadySimulator, x_dof::AbstractArray)::Tuple{Real, Real}
    err = Vector{Float64}()
    @inbounds for (_, pipe) in ref(ss, :pipe)
        eqn_no = pipe["dof"]
        f = x_dof[eqn_no]
        fr_node = pipe["fr_node"]  
        to_node = pipe["to_node"]
        fr_dof = ref(ss, :node, fr_node, "dof")
        to_dof = ref(ss, :node, to_node, "dof")

        R1 = nominal_values(ss, :mach_num)^2 / (nominal_values(ss, :euler_num) * pipe["area"]^2)  
        R2 = nominal_values(ss, :mach_num)^2 / ( nominal_values(ss, :euler_num) * nominal_values(ss, :froude_num)^2 )
        beta = pipe["friction_factor"]  / (2 * pipe["diameter"])

        R1_bar = R1 / nominal_values(ss, :euler_num)
        R2_bar = R2 * nominal_values(ss, :euler_num)
        # 8 degree inclination = sin(theta) approx 0.14
        # 2 degree = sin(theta) approx 0.034
        sin_incline = 0.065
        p_fr = dof_to_pressure(x_dof[fr_dof])
        p_to = dof_to_pressure(x_dof[to_dof])
        
        var1 = ((beta * R1_bar) / (sin_incline * R2_bar)) * f * abs(f)
        var2 = (R1_bar) * (f^2)
        var = abs((var2 - var1) * log( abs(p_fr^2 - var1) / abs(p_to^2 - var1) ) 
        - (2) * var2 * log( abs(p_fr) / abs(p_to) ) - 2 * pipe["length"] * beta * R1_bar * f * abs(f))
        push!(err, var)
    end
    err_max  = maximum(err)
    err_rms = sqrt(sum(err.^2)/ length(err))
    return err_max, err_rms
end
