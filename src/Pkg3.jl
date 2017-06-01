module Pkg3

using TOML

# Pkg3 package loaders

abstract type Loader end

struct LoadInstalled <: Loader
    depot::String
end

struct AutoInstall <: Loader
    depot::String
end

function Base.load_hook(loader::LoadInstalled, name::String, found)
    # load and parse the package manifest
    manifest = TOML.parsefile("Manifest.toml")
    # pull the package UUID and verion hash out of the manifest
    packages = haskey(manifest, "package")  ? manifest["package"]  : return nothing
    package  = haskey(packages, name)       ? packages[name]       : return nothing
    uuid     = haskey(package, "uuid")      ? package["uuid"]      : return nothing
    hash     = haskey(package, "hash-sha1") ? package["hash-sha1"] : return nothing
    # see if uuid/hash is installed in the current depot:
    path = joinpath(loader.depot, "packages", uuid, hash, "src", "$name.jl")
    return isfile(path) ? path : nothing
end

function Base.load_hook(loader::AutoInstall, name::String, found)
    return found
end

user_depot() = abspath(homedir(), ".julia")

push!(LOAD_PATH, LoadInstalled(user_depot()))
unshift!(LOAD_PATH, AutoInstall(user_depot()))

end # module
