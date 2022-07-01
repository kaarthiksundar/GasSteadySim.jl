function run_simulator!(ss::SteadySimulator; 
    method::Symbol=:newton,
    iteration_limit::Int64=2000, 
    kwargs...)
    

    ## Stage 1
    A_ns, A_alpha_ns, A_slack, beta_vec, val_ns, pi_slack_vec = assemble_mat(ss)

    # println(A_ns, "\n", A_alpha_ns, "\n", A_slack, "\n", beta_vec, "\n", val_ns)

    num_ns =  ss.new_ref[:total_nonslack_vertices]
    num_edges = ss.new_ref[:total_edges]

    j_vec, i_vec, k_alpha_vec = findnz(A_alpha_ns)
    A_alpha_ns_transpose = sparse(i_vec, j_vec, k_alpha_vec, num_edges, num_ns)

    # @show A_alpha_ns, A_alpha_ns_transpose, pi_slack_vec, beta_vec, val_ns

    # return A_ns, A_alpha_ns, num_ns, num_edges

    ## Stage 2
    gamma_max = 1e7  #1e4 gaslib40

    # while true
    #     status, W1 , W2 = find_gamma_W(num_ns, num_edges, A_ns, A_alpha_ns, beta_vec, gamma_max)
    #     if status == MathOptInterface.OPTIMAL
    #         println("Success ", gamma_max, "\n")
    #         gamma_max = 10 * gamma_max 
    #         continue
    #     else
    #         println("Failed for ", gamma_max, " ", status, "\n")
    #         gamma_max = gamma_max/10
    #         break
    #     end
    # end 

    println(gamma_max)
    status, W1, W2 = find_gamma_W(num_ns, num_edges, A_ns, A_alpha_ns, beta_vec, gamma_max)
    # note that premultiplying matrix for monotone operator is not just W1, W2  but W*D , D is block diagonal with I, beta_vec 
    W2_final =  -W2 * Diagonal(1.0 ./ beta_vec)

    # @show W2_final

 

    ## Stage 3
    ## next  block required to compute residual
    residual_nonslack = zeros(Float64, num_ns) 
    residual_edges = zeros(Float64, num_edges)

    #here do VI soln algo in loop
    x_old =  -ones(Float64, num_ns + num_edges)
    x_int = zeros(Float64, num_ns + num_edges)
    x_new = zeros(Float64, num_ns + num_edges)
    residual_full = zeros(Float64, num_ns + num_edges)

    #Stage 4
    # Estimate Lipschitz const L, verify operator is monotone

#     for i = 1:10
#         y1 =  rand(Float64, num_ns + num_edges) .+ 3.0
#         y2 =  rand(Float64, num_ns + num_edges)
# # 
#         assemble_residual!(ss, y1, residual_edges, residual_nonslack, A_ns, A_alpha_ns_transpose, A_slack, val_ns, beta_vec, pi_slack_vec, W1, W2_final)
#         x_old .= [residual_nonslack; residual_edges]
#         assemble_residual!(ss, y2, residual_edges, residual_nonslack, A_ns, A_alpha_ns_transpose, A_slack, val_ns, beta_vec, pi_slack_vec, W1, W2_final)
#         x_new  .= [residual_nonslack; residual_edges]
#         L = norm(x_old - x_new) / norm(y1 - y2)
#         m = dot(x_old - x_new, y1 - y2)
#         println("L: ", L, " m: ", m)
#     end
    


    ## Stage 5
    ## Find zero of monotone operator by extra gradient method
    tau = 1e-4 # tau < 1/L  #1e-4, 600 iters was enough for 8-node. 1e-6, 80e3 for gaslib40
    lb = 3e-1
    crossings = zeros(Int64, num_edges)

    for i = 1:100
        assemble_residual!(ss, x_old, residual_edges, residual_nonslack, A_ns, A_alpha_ns_transpose, A_slack, val_ns, beta_vec, pi_slack_vec, W1, W2_final)
        # println("ns: ", residual_nonslack, "\n", "edges: ", residual_edges)
        residual_full .= [residual_nonslack; residual_edges]
        # 

        x_int .= x_old - tau * [residual_nonslack; residual_edges] 
        # @show x_int
        # _project_onto_domain_1!(ss, x_int, x_old, lb, gamma_max)
         _project_onto_domain_2!(ss, x_int, x_old, crossings, lb, gamma_max)
        # @show x_int
        #
        assemble_residual!(ss, x_int, residual_edges, residual_nonslack, A_ns, A_alpha_ns_transpose, A_slack, val_ns, beta_vec, pi_slack_vec, W1, W2_final)
        x_new .= x_old .-  tau * [residual_nonslack; residual_edges]
        # @show x_new
        _project_onto_domain_2!(ss, x_new, x_old, crossings, lb, gamma_max)
        detect_repeated_crossings(crossings)
        # @show x_new

        if i%20 == 0
            # println("ns: ", residual_nonslack, "\n", "edges: ", residual_edges)
            @show crossings
            println(norm(residual_full), "\n")
            println("xdiff: ", norm(x_old - x_new), "\n")
        end
        
        println(norm(x_new - x_old), "\n")
        x_old .= x_new
        # @show i, x_old[num_ns+1]
    end
    # @show x_old

    # Now map x to potential, edge flow in ref (not new_ref)
    # post process to get compressor edge flows, delivery vertex pressures
    # here operator is Lipschitz over a finite domain, not all of Rn. Does this mean extra gradient method doesn't apply though it works ?
    # can we use extragradient mthd with projection to solve VI instead of iteration to find 0 ? Then we can show uniqueness of flow solution.


end


function find_gamma_W(num_ns::Int64, num_edges::Int64, A_ns::SparseArrays.SparseMatrixCSC, A_alpha_ns::SparseArrays.SparseMatrixCSC, beta_vec::AbstractArray, gamma::Real)
    xmax = 1.0
    model = Model(SCS.Optimizer)
    set_silent(model)
    @variable(model, Wa[i= 1:num_ns, j = 1:num_ns], PSD)
    @variable(model, Wb[i= 1:num_edges, j = 1:num_edges], PSD)
    @variable(model, 1e-1 <= x <= xmax)
    @constraint(model, (Wa .+ x) * A_ns .- A_alpha_ns * Diagonal(1.0 ./ beta_vec) * (Wb .+ x)   .== zeros(num_ns, num_edges) )
    @constraint(model, [i=1:num_edges], (gamma + 1) * Wb[i, i] + x >=  gamma * sum(Wb[i, j] for j=1:num_edges) )
    # bounds on norms of Wa, Wb lead to infeasibility, so no norm constraint
    @objective(model, Max, x)
    optimize!(model)

    @show termination_status(model)
    if termination_status(model) != OPTIMAL
        @warn("Failed")
        return termination_status(model), nothing, nothing
    end
    @show value(x)
    W1 = value(Wa) .+ value(x)
    W2 = value(Wb) .+ value(x)
    println(eigmin(W1), " ", eigmin(W2), "\n")

    return termination_status(model), W1, W2

end

function _project_onto_domain_1!(ss::SteadySimulator, x_int::AbstractArray, x_old::AbstractArray, epsilon::Real, gamma::Real)
    num_ns =  ss.new_ref[:total_nonslack_vertices]
    num_edges = ss.new_ref[:total_edges]
    # sign_vec = sign.(residual_full[num_ns+1:num_ns+num_edges])
    # residual_full[num_ns+1:num_ns+num_edges] .= max.( epsilon, min.( abs.(residual_full[num_ns+1:num_ns+num_edges]), gamma * epsilon ) ) .* sign_vec
    

    # crossings flag for each edge, all 0
    for i = 1 : num_edges
        if x_int[i + num_ns] > gamma * epsilon
            x_int[i + num_ns] = gamma * epsilon
        elseif x_int[i + num_ns] < -gamma * epsilon
            x_int[i + num_ns] = -gamma * epsilon
        elseif x_int[i + num_ns] < epsilon < x_old[i + num_ns]
            x_int[i + num_ns] = epsilon
        elseif x_old[i + num_ns] < -epsilon < x_int[i + num_ns]
            x_int[i + num_ns] = -epsilon
        end
    end

    return

end

function _project_onto_domain_2!(ss::SteadySimulator, x_new::AbstractArray, x_old::AbstractArray, crossings::AbstractArray, epsilon::Real, gamma::Real)
    num_ns =  ss.new_ref[:total_nonslack_vertices]
    num_edges = ss.new_ref[:total_edges]
    # sign_vec = sign.(residual_full[num_ns+1:num_ns+num_edges])
    # residual_full[num_ns+1:num_ns+num_edges] .= max.( epsilon, min.( abs.(residual_full[num_ns+1:num_ns+num_edges]), gamma * epsilon ) ) .* sign_vec
    

    # crossings flag for each edge, all 0
    for i = 1 : num_edges
        if x_new[i + num_ns] > gamma * epsilon
            x_new[i + num_ns] = gamma * epsilon

        elseif x_new[i + num_ns] < -gamma * epsilon
            x_new[i + num_ns] = -gamma * epsilon

        elseif (-epsilon < x_new[i + num_ns] < 0 )  && (epsilon < x_old[i + num_ns] < gamma*epsilon)
            crossings[i] += 1
            x_new[i + num_ns] = -epsilon

        elseif (0 < x_new[i + num_ns] < epsilon )  && (epsilon < x_old[i + num_ns] < gamma*epsilon)
            crossings[i] += 1
            x_new[i + num_ns] = epsilon

        elseif (-epsilon < x_new[i + num_ns] < 0 )  && (-gamma*epsilon < x_old[i + num_ns] < -epsilon)
            crossings[i] += 1
            x_new[i + num_ns] = -epsilon

        elseif (0 < x_new[i + num_ns] < epsilon )  && (-gamma*epsilon < x_old[i + num_ns] < -epsilon)
            crossings[i] += 1
            x_new[i + num_ns] = epsilon

        end
    end

    return

end


function detect_repeated_crossings(crossings::AbstractArray)

    if sum( [crossings == 2 for i = 1:length(crossings)] ) > 0
        @error "repeated crossing between domains. Solution may not exist"
    end
    return
end

function _project_onto_domain!(ss::SteadySimulator, x_new::AbstractArray, epsilon::Real, gamma::Real)
    num_ns =  ss.new_ref[:total_nonslack_vertices]
    num_edges = ss.new_ref[:total_edges]
    # sign_vec = sign.(residual_full[num_ns+1:num_ns+num_edges])
    # residual_full[num_ns+1:num_ns+num_edges] .= max.( epsilon, min.( abs.(residual_full[num_ns+1:num_ns+num_edges]), gamma * epsilon ) ) .* sign_vec
    

    # crossings flag for each edge, all 0
    for i = num_ns + 1 : num_ns + num_edges
        if x_new[i] > gamma * epsilon
            x_new[i] = gamma * epsilon
        elseif x_new[i] < -gamma * epsilon
            x_new[i] = -gamma * epsilon
        elseif 0 < x_new[i] < epsilon
            crossings[i] += 1
            x_new[i] = epsilon
        elseif -epsilon < x_new[i] < 0
            crossings[i] += 1
            x_new[i] = -epsilon
        end
    end

    return

end
# function _create_initial_guess_dof!(ss::SteadySimulator)::Array
#     ndofs = length(ref(ss, :dof))
#     x_guess = 0.5 * ones(Float64, ndofs) 
#     dofs_updated = 0

#     components = [:node, :pipe, :compressor, 
#         :control_valve, :valve, 
#         :resistor, :loss_resistor, :short_pipe]

#     for component in components 
#         for (i, val) in get(ss.initial_guess, component, [])
#             x_guess[ref(ss, component, i, "dof")] = val 
#             dofs_updated += 1
#         end 
#     end 
#     return x_guess
# end