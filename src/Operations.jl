module Operations

using ..VersionSpecs
using Base.Random: UUID
using TOML

depots() = Base.Loading.DEPOTS

function registries(depot::String)
	d = joinpath(depot, "registries")
	filter!(readdir(d)) do r
		isfile(joinpath(d, r, "registry.toml")) &&
		isfile(joinpath(d, r, "packages.toml"))
	end
end

function parse_packages(depot::String, registry::String)
	packages_toml = joinpath(depot, "registries", registry, "packages.toml")
	convert(Dict{UUID,Dict{Symbol,String}}, TOML.parsefile(packages_toml))
end

function find_registered(names::Vector{String})
	# name --> uuid --> paths
	paths = Dict{String,Dict{UUID,Vector{String}}}()
	for depot in depots(), registry in registries(depot)
		for (uuid, info) in parse_packages(depot, registry)
			name = info[:name]
			name in names || continue
			if !haskey(paths, name)
				paths[name] = Dict{UUID,Vector{String}}()
			end
			if !haskey(paths[name], uuid)
				paths[name][uuid] = String[]
			end
			path = abspath(registry, info[:path])
			push!(paths[name][uuid], path)
		end
	end
	return paths
end

function add(names::Dict{String,VersionSpec})

end

end # module
