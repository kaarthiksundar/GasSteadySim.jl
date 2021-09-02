
function _build_ig(data::Dict{String,Any})::Dict{Symbol,Any}
    ig = Dict{Symbol,Any}()

    ig[:node] = Dict() 
    ig[:pipe] = Dict()
    ig[:compressor] = Dict()

    for (i, value) in get(data, "nodal_pressure", [])
        id = parse(Int64, i)
        ig[:node][id] = value
    end 

    for (i, value) in get(data, "pipe_flow", [])
        id = parse(Int64, i) 
        ig[:pipe][id] = value 
    end 

    for (i, value) in get(data, "compressor_flow", [])
        id = parse(Int64, i)
        ig[:compressor][id] = value
    end 

    (isempty(ig[:node])) && (@info "initial guess for nodal pressure not provided") 
    (isempty(ig[:pipe])) && (@info "initial guess for pipe flow not provided")
    (isempty(ig[:compressor])) && (@info "initial guess for compressor flow not provided")

    return ig
end 