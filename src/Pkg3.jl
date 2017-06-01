module Pkg3

struct Loader
	env::Dict{String}
	depot::String
end

function Base.load_hook(prefix::Loader, name::String, found)
	
end

end # module
