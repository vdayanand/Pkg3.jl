module Pkg3

include("Types.jl")
include("Operations.jl")
include("Loading.jl")

import .Types: VersionSpec, @vs_str
export VersionSpec, @vs_str

import .Loading: Loader, LoadInstalled
export Loader, LoadInstalled

end # module
