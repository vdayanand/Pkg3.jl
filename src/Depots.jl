module Depots

using TOML

function registries(depot::String)
	d = joinpath(depot, "registries")
	filter!(readdir(d)) do r
		isfile(joinpath(d, r, "registry.toml")) &&
		isfile(joinpath(d, r, "packages.toml"))
	end
end

function parse_packages(depot::String, registry::String)
	packages_toml = joinpath(depot, "registries", registry, "packages.toml")
	packages = convert(Dict{String,String}, TOML.parsefile(packages_toml))
end

function package_name(depot::String, registry::String, uuid::String)
	packages = parse_packages(depot, registry)
	get(packages, uuid, nothing)
end

function package_uuid(depot::String, registry::String, name::String)
	packages = Dict(n => u for (u, n) in parse_packages(depot, registry))
	get(packages, name, nothing)
end

function install(depot::String, name::String, version::VersionNumber)
	for registry in registries(depot)
		uuid = 
	end
end

end # module
