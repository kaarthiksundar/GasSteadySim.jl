function _si_to_english!(data::Dict{String,Any},
    params::Dict{Symbol,Any}, nominal_values::Dict{Symbol,Any})

    rescale_mass_flow = x -> x / get_mmscfd_to_kgps_conversion_factor(params)
    rescale_mass_flux = x -> x / get_mmscfd_to_kgps_conversion_factor(params)
    rescale_pressure = x -> pascal_to_psi(x)
    rescale_length = x -> m_to_miles(x)
    rescale_density = x -> x
    rescale_diameter = x -> m_to_inches(x)
    rescale_area = x -> sq_m_to_sq_inches(x)
    function rescale_compressor_boundary_conditions!(type, value)
        (type == 1) && (value = rescale_pressure(value))
        (type == 2) && (value = rescale_mass_flow(value))
    end 

    rescale_functions = [rescale_mass_flow, rescale_mass_flux, 
        rescale_time, rescale_pressure, rescale_length, rescale_density, 
        rescale_diameter, rescale_area, rescale_compressor_boundary_conditions!]
    
    _rescale_data!(data, params, rescale_functions)
end 