module Pkg3

include("Loading.jl")
include("Depots.jl")
include("Environments.jl")

import .Loading: Loader, LoadInstalled

export Loader, LoadInstalled

end # module
