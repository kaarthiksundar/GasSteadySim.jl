
function _initialize_solution(data::Dict{String,Any})::Dict{String,Any}
    sol = Dict{String,Any}()
    sol["nodal_pressure"] = Dict{Int64,Float64}()
    sol["nodal_density"] = Dict{Int64,Float64}()
    sol["pipe_flow"] = Dict{Int64,Float64}()
    (haskey(data, "compressors")) && (sol["compressor_flow"] = Dict{Int64,Float64}())
    (haskey(data, "control_valves")) && (sol["control_valve_flow"] = Dict{Int64,Float64}())
    (haskey(data, "valves")) && (sol["valve_flow"] = Dict{Int64,Float64}())
    (haskey(data, "resistors")) && (sol["resistor_flow"] = Dict{Int64,Float64}())
    (haskey(data, "loss_resistors")) && (sol["loss_resistor_flow"] = Dict{Int64,Float64}())
    (haskey(data, "short_pipes")) && (sol["short_pipe_flow"] = Dict{Int64,Float64}())
    return sol
end