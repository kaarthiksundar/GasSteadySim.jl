function populate_lp_model!(ss::SteadySimulator)
    model = ss.lp_relax
    m = model.model
    var = model.variables 
    sol = model.sol 
    residual = model.residuals
    
end 