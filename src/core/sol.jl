
function _initialize_solution(data::Dict{String,Any})::Dict{String,Any}
    sol = Dict{String,Any}()
    sol["nodal_pressure"] = Dict{Int64,Float64}()
    sol["pipe_flow"] = Dict{Int64,Float64}()
    if haskey(data, "compressors")
        sol["compressor_flow"] = Dict{Int64,Float64}()
    end
    if haskey(data, "control_valves") 
        sol["control_valve_flow"] = Dict{Int64,Float64}()
    end 
    return sol
end