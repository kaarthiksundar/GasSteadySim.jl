module GasSteadySim

import JSON
using NLsolve
using SparseArrays
using LineSearches
# using LinearAlgebra # not needed if not doing fixed point

include("io/json.jl")
include("io/data_utils.jl")

include("unit_conversion/unit_convertor_utils.jl")
include("unit_conversion/to_si.jl")
include("unit_conversion/to_english.jl")
include("unit_conversion/to_pu.jl")
include("unit_conversion/unit_convertors.jl")

include("core/eos.jl")
include("core/types.jl")
include("core/ref.jl")
include("core/ig.jl")
include("core/bc.jl")
include("core/sol.jl")
include("core/bounds.jl")
include("core/initialize_ss.jl")
include("core/assemble.jl")
include("core/run_ss.jl")
include("core/output.jl")
include("io/writer.jl")
include("core/export.jl")

end # module
