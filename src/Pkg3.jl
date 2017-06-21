module Pkg3

include("Types.jl")
include("Operations.jl")
include("Loading.jl")

import .Types: VersionSpec, @vs_str
import .Operations: add
import .Loading: Loader, LoadInstalled

end # module
