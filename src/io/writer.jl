function write_output(ss::SteadySimulator;
    output_path::AbstractString="./", 
    output_file::AbstractString="output.json" 
    )

    solution = ss.sol 
    output_string = output_path * output_file

    open(output_string, "w") do f 
        JSON.print(f, solution, 2)
    end 
end 