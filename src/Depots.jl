module Depots

using TOML

function registries(depot::String)
	d = joinpath(depot, "registries")
	filter!(readdir(d)) do r
		isfile(joinpath(d, r, "registry.toml")) &&
		isfile(joinpath(d, r, "packages.toml"))
	end
end

function _parse_packages(depot::String, registry::String)
	packages_toml = joinpath(depot, "registries", registry, "packages.toml")
	packages = convert(Dict{String,String}, TOML.parsefile(packages_toml))
end

function package_name(depot::String, registry::String, uuid::String)
	packages = _parse_packages(depot, registry)
	packages[uuid] # throws KeyError => PackageNotFound(uuid)
end

function package_uuid(depot::String, registry::String, name::String)
	packages = Dict(n => u for (u, n) in _parse_packages(depot, registry))
	packages[name] # throws KeyError => PackageNotFound(name)
end

function install(depot::String, name::String, version::VersionNumber)

end

end # module
