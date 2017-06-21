module Pkg3

include("Types.jl")
include("Operations.jl")
include("Loading.jl")

import .Types: VersionSpec, @vs_str
import .Operations: add, find_config, find_manifest, user_depot, depots, registries
import .Loading: Loader, LoadInstalled

export @vs_str

end # module
