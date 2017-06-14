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

function install(pkg::String, ver::VersionNumber)
	uuids = String[]
	for depot in Base.Loading.DEPOTS, registry in registries(depot)
		for (uuid, name) in parse_packages(depot, registry)
			name == pkg && push!(uuids, uuid)
		end
	end
	return uuids
end

end # module
