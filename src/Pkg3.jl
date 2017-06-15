module Pkg3

include("Loading.jl")
include("Operations.jl")

import .Loading: Loader, LoadInstalled

export Loader, LoadInstalled

end # module
