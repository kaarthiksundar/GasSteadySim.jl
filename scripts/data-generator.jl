using JSON
using GasModels
using Base.Filesystem

get_simulation_params(; T::Float64 = 288.70599999999996) = 
    Dict{String,Any}(
    "Temperature (K):" => T,
    "Gas specific gravity (G):" => 0.6,
    "Specific heat capacity ratio" => 1.4,
    "units (SI = 0, standard = 1)" => 0.0
)


function generate_steady_sim_inputs_from_gaslib(file::AbstractString, output_folder::AbstractString)

    data = parse_file("GasLib-benchmarks/" * file)

    (!isdir(output_folder)) && (mkdir(output_folder)) 

    network_data = Dict{String,Any}(
        "nodes" => Dict{String,Any}(), 
        "pipes" => Dict{String,Any}(), 
        "compressors" => Dict{String,Any}()
    )

    bc = Dict{String,Any}(
        "boundary_compressor" => Dict{String,Any}(),
        "boundary_nonslack_flow" => Dict{String,Any}(),
        "boundary_pslack" => Dict{String,Any}()
    )

    withdrawal = Dict{String,Any}()

    for (i, _) in data["junction"]
        withdrawal[i] = 0.0
    end 

    for (_, receipt) in data["receipt"]
        node = string(receipt["junction_id"])
        withdrawal[node] -= receipt["injection_nominal"] * data["base_flow"]
    end 

    for (_, delivery) in data["delivery"]
        node = string(delivery["junction_id"])
        withdrawal[node] += delivery["withdrawal_nominal"] * data["base_flow"]
    end 

    withdrawal_values = [v for (k, v) in withdrawal]
    max_injection = minimum(withdrawal_values)
    nodes_with_max_injection = sort([k for (k,v) in withdrawal if v == max_injection])
    slack_node = nodes_with_max_injection[1]

    for (i, node) in data["junction"]
        network_data["nodes"][i] = Dict{String,Any}(
            "node_id" => parse(Int, i),
            "node_name" => "n" * i,
            "x_coord" => node["lat"],
            "y_coord" => node["lon"], 
            "min_pressure" => node["p_min"] * data["base_pressure"], 
            "max_pressure" => node["p_max"] * data["base_pressure"]
        )
        if i == string(slack_node)
            network_data["nodes"][i]["slack_bool"] = 1
            bc["boundary_pslack"][i] = 5000000
        else 
            network_data["nodes"][i]["slack_bool"] = 0
            bc["boundary_nonslack_flow"][i] = withdrawal[i]
        end
    end

    for (i, pipe) in data["pipe"]
        network_data["pipes"][i] = Dict{String,Any}(
            "pipe_id" => parse(Int, i),
            "pipe_name" => "p" * i,
            "from_node" => pipe["fr_junction"],
            "to_node" => pipe["to_junction"],
            "diameter" => pipe["diameter"] * data["base_diameter"],
            "length" => pipe["length"] * data["base_length"],
            "friction_factor" => pipe["friction_factor"],
        )
    end 

    for (i, compressor) in data["compressor"]
        network_data["compressors"][i] = Dict{String,Any}(
            "comp_id" => parse(Int, i),
            "comp_name" => "c" * i,
            "from_node" => compressor["fr_junction"],
            "to_node" => compressor["to_junction"]
        )
        bc["boundary_compressor"][i] = Dict{String,Any}(
            "control_type" => 0, "value" => 1.5
        )
    end 

    open(output_folder * "params.json", "w") do f 
        JSON.print(f, Dict("simulation_params" => get_simulation_params()), 2)
    end
 
    open(output_folder * "network.json", "w") do f 
        JSON.print(f, network_data, 2)
    end

    open(output_folder * "bc.json", "w") do f 
        JSON.print(f, bc, 2)
    end

end 

generate_steady_sim_inputs_from_gaslib("GasLib-40.zip", "GasLib-40/")
generate_steady_sim_inputs_from_gaslib("GasLib-135.zip", "GasLib-135/")