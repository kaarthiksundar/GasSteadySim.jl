
function initialize_solution(data::Dict{String,Any})::Dict{String,Any}
    sol = Dict{String,Any}()
    sol["nodal_pressure"] = Dict{Int64,Float64}()
    sol["pipe_flow"] = Dict{Int64,Float64}()
    sol["compressor_flow"] = Dict{Int64,Float64}()
    
    return sol
end