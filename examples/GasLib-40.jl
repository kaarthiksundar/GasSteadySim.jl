using GasSteadySim
using Distributions
using DelimitedFiles

function randomize_compressor_ratio!(ss::SteadySimulator, lb::Float64, ub::Float64)
	num_compressors = length(ref(ss, :compressor))

	var = rand(Uniform(lb, ub), num_compressors)
	for i = 1:num_compressors
		ss.boundary_conditions[:compressor][i]["val"] = var[i]
		ss.ref[:compressor][i]["c_ratio"] = var[i]
	end

end

function randomize_nonslack_injections!(ss::SteadySimulator, lb::Float64, ub::Float64)

	s = Vector{Int32}[]
	for (i, _) in ref(ss, :node)
		if ref(ss, :node, i, "is_slack") == 1
			push!(s, i)
		end
	end

	num_nodes = length(ref(ss, :node))
	var = rand(Uniform(lb, ub), num_nodes)
	for i = 1:num_nodes
		if i in s
			continue
		end
		q = ss.boundary_conditions[:node][i]["val"]
		ss.boundary_conditions[:node][i]["val"] = var[i] * q
		ss.ref[:node][i]["withdrawal"] = var[i] * q
	end

end

function write_bc_fail(ss::SteadySimulator, output_file::AbstractString;
    output_path::AbstractString="./")

	data = Dict{String, Any}()
	data["boundary_pslack"] = Dict{Int64, Float64}()
	data["boundary_compressor"] = Dict{Int64, Any}()
	data["boundary_nonslack_flow"] = Dict{Int64, Float64}()

	for (i, _) in ref(ss, :node)
		if ref(ss, :node, i, "is_slack") == 1
			data["boundary_pslack"][i] = ref(ss, :node, i, "pressure") * nominal_values(ss, :pressure)
		else
			data["boundary_nonslack_flow"][i] = ref(ss, :node, i, "withdrawal") * nominal_values(ss, :mass_flow) 
		end
	end


	num_compressors = length(ref(ss, :compressor))
	for i = 1:num_compressors
		data["boundary_compressor"][i] = Dict{String, Any}()
		data["boundary_compressor"][i]["control_type"] = 0
		data["boundary_compressor"][i]["value"] = ref(ss, :compressor, i, "c_ratio") 
	end


    output_string = output_path * output_file

    open(output_string, "w") do f 
        JSON.print(f, data, 2)
    end 
end 


file = "examples/data/GasLib-40/"
ss = initialize_simulator(file, eos=:simple_cnga, initial_guess_filename="")
solver_return = run_simulator!(ss; jacobian_type = :in_place)

println(solver_return.status)
