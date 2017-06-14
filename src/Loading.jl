module Loading

using TOML

# Pkg3 package loaders

abstract type Loader end

struct LoadInstalled <: Loader
    depot::String
end

function Base.load_hook(loader::LoadInstalled, name::String, found)
    # load and parse the package manifest
    manifest = TOML.parsefile("Manifest.toml")
    # pull the package UUID and verion hash out of the manifest
    packages = haskey(manifest, "package")  ? manifest["package"]  : return found
    package  = haskey(packages, name)       ? packages[name]       : return found
    uuid     = haskey(package, "uuid")      ? package["uuid"]      : return found
    hash     = haskey(package, "hash-sha1") ? package["hash-sha1"] : return found
    # relative path for package and version
    path = joinpath(uuid, hash, "src", "$name.jl")
    # see if previously found matches expected suffix
    found isa String && endswith(found, path) && return found
    # see if package and version is installed in the current depot
    path = joinpath(loader.depot, "packages", path)
    # if it is, return it, otherwise pass through
    return isfile(path) ? path : found
end

user_depot() = abspath(homedir(), ".julia")

push!(LOAD_PATH, LoadInstalled(user_depot()))

@eval Base module Loading; DEPOTS = []; end

push!(Base.Loading.DEPOTS, user_depot())

end # module
