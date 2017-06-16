module Pkg3

include("Loading.jl")
include("VersionSpecs.jl")
include("Operations.jl")

import .Loading: Loader, LoadInstalled
export Loader, LoadInstalled

import .VersionSpecs: VersionSpec, @vs_str
export VersionSpec, @vs_str

end # module
