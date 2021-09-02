function _si_to_pu!(data::Dict{String,Any},
    params::Dict{Symbol,Any}, nominal_values::Dict{Symbol,Any})

    rescale_mass_flow = x -> x/nominal_values[:mass_flow]
    rescale_mass_flux = x -> x/nominal_values[:mass_flux]
    rescale_pressure = x -> x/nominal_values[:pressure]
    rescale_length = x -> x/nominal_values[:length]
    rescale_density = x -> x/nominal_values[:density]
    rescale_diameter = x -> x/nominal_values[:length]
    rescale_area = x -> x/nominal_values[:area]


    function rescale_compressor_boundary_conditions!(value)
        
        type = value["control_type"]
        val =  value["value"]
        (type == 1) && (value["value"] = rescale_pressure(val))
        (type == 2) && (value["value"] = rescale_mass_flow(val))
        
    end  

    rescale_functions = [rescale_mass_flow, rescale_mass_flux, rescale_pressure, rescale_length, rescale_density, 
        rescale_diameter, rescale_area, rescale_compressor_boundary_conditions!]
    
    _rescale_data!(data, params, rescale_functions)
end 