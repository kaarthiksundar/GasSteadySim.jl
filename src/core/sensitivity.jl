
function initialize_for_sensitivity_computations(ss::SteadySimulator)::AbstractArray

	ndofs = length(ref(ss, :dof))
    RHS = zeros(ndofs, ref(ss, :total_control_vars))
    return RHS 

end

function assemble_rhs_corresponding_to_single_control_var!(ss::SteadySimulator,  i::Int64,  RHS::AbstractArray, x_dof::AbstractArray)

	component, index = ref(ss, :control_vars, i)
	if component == :node

		RHS[ref(ss, :node, index, "dof"), i] = 1.0
		
	elseif component == :compressor

        fr_node = ref(ss, :compressor, index, "fr_node")

        RHS[ref(ss, :compressor, index, "dof"), i] = -1 * x_dof[ref(ss, :node, fr_node, "dof")]

	end
	return
end

function assemble_rhs_with_all_control_var!(ss::SteadySimulator,  RHS::AbstractArray, x_dof::AbstractArray)
	for i = 1 : ref(ss, :total_control_vars)
        assemble_rhs_corresponding_to_single_control_var!(ss, i,  RHS, x_dof)
    end
    return
end


function calculate_sensitivities!(ss::SteadySimulator, J::AbstractArray, RHS::AbstractArray)

	lufact = lu(J)
	for i = 1 : ref(ss, :total_control_vars)
		x = lufact \ RHS[:, i]
		RHS[:, i] .= x
	end
	return
end


function convert_sensitivities_to_given_units(ss::SteadySimulator, RHS::AbstractArray)::AbstractArray

	ndofs = length(ref(ss, :dof))
	Sensitivity_matrix = zeros(ndofs, ref(ss, :total_control_vars))

	for i = 1: length(ref(ss, :dof))
		num_scaling = 0.0
		denom_scaling = 0.0
		sym, local_id = ref(ss, :dof, i)
		if sym == :node
			num_scaling = (ref(ss, :is_pressure_node, local_id)) ? ss.nominal_values[:pressure] : get_potential(ss, ss.nominal_values[:pressure])

		elseif sym == :pipe
			num_scaling = ss.nominal_values[:mass_flow]
		elseif sym == :compressor
			num_scaling = 1.0
		else
			@error "Not implemented anything other than nodal pressures/injections, pipe flows and compressor ratios"
		end

		for col = 1 : ref(ss, :total_control_vars)
			component, index = ref(ss, :control_vars, col)

			if component == :node 
				if ref(ss, :node, index, "is_slack") == 1
					denom_scaling = (ref(ss, :is_pressure_node, local_id)) ? ss.nominal_values[:pressure] : get_potential(ss, ss.nominal_values[:pressure])
				else
					denom_scaling = ss.nominal_values[:mass_flow]
				end 

			elseif component == :compressor
				denom_scaling = 1.0
			end

			Sensitivity_matrix[i, col] = (num_scaling/denom_scaling) * RHS[i, col]

		end
	end

	return Sensitivity_matrix

end