module Operations

using TOML
using Base.Random: UUID
using Base: LibGit2
using Pkg3.Types

user_depot() = abspath(homedir(), ".julia")
depots() = Base.Loading.DEPOTS

function registries(depot::String)
    d = joinpath(depot, "registries")
    regs = filter!(readdir(d)) do r
        isfile(joinpath(d, r, "registry.toml")) &&
        isfile(joinpath(d, r, "packages.toml"))
    end
    return map(reg->joinpath(depot, "registries", reg), regs)
end
registries() = [r for d in depots() for r in registries(d)]

function parse_packages(registry::String)
    packages_toml = joinpath(registry, "packages.toml")
    convert(Dict{UUID,Dict{Symbol,String}}, TOML.parsefile(packages_toml))
end

function find_registered(names::Vector{String})
    # name --> uuid --> paths
    where = Dict{String,Dict{UUID,Vector{String}}}()
    for registry in registries()
        for (uuid, info) in parse_packages(registry)
            name = info[:name]
            name in names || continue
            if !haskey(where, name)
                where[name] = Dict{UUID,Vector{String}}()
            end
            if !haskey(where[name], uuid)
                where[name][uuid] = String[]
            end
            path = abspath(registry, info[:path])
            push!(where[name][uuid], path)
        end
    end
    return where
end

include("libgit2_discover.jl")

default_env() = get(ENV, "JULIA_ENV", nothing)
find_env() = find_env(default_env())

function find_env(env::String)
    names = ["JuliaConfig.toml", "Config.toml"]
    if isempty(env) || env == "."
        # find project environment...
        path = LibGit2.discover(ceiling = homedir())
        repo = LibGit2.GitRepo(path)
        work = LibGit2.workdir(repo)
        for name in names
            path = abspath(work, name)
            isfile(path) && return path
        end
        return abspath(work, names[end])
    elseif startswith(env, "/") || startswith(env, "./")
        # path to config file or project directory
        splitext(env)[2] == ".toml" && return abspath(env)
        for name in names
            path = abspath(env, name)
            isfile(path) && return path
        end
        return abspath(env, names[end])
    else # named environment
        for depot in depots()
            path = joinpath(depot, "environments", env, "Config.toml")
            isfile(path) && return path
        end
        return joinpath(user_depot(), "environments", env, "Config.toml")
    end
end

function find_env(::Void)
    # look for a default environment
    defaults = [
        "v$(VERSION.major).$(VERSION.minor).$(VERSION.patch)",
        "v$(VERSION.major).$(VERSION.minor)",
        "v$(VERSION.major)",
        "default",
    ]
    for depot in depots(), env in defaults
        path = joinpath(depot, "environments", env, "Config.toml")
        isfile(path) && return path
    end
    env = VERSION.major == 0 ? defaults[2] : defaults[3]
    return joinpath(user_depot(), "environments", env, "Config.toml")
end

function find_manifest(env::Union{Void,String}=default_env())
    config = find_env(env)
    isfile(config) && return abspath(TOML.parsefile(config)["manifest"])
    names = ["JuliaManifest.toml", "Manifest.toml"]
    dir = dirname(config)
    for name in names
        path = joinpath(dir, name)
        isfile(path) && return path
    end
    return joinpath(dir, names[end])
end

function add(pkgs::Dict{String})
    names = sort!(collect(keys(pkgs)))
    where = find_registered(names)
    # check for ambiguous package names
    uuids = let ambig = false
        for name in names
            length(where[name]) == 1 && continue
            msg = "$name is ambiguous, it could refer to:\n"
            for (i, (uuid, paths)) in enumerate(sort!(collect(where[name]), by=first))
                msg *= " [$i] $uuid"
                for path in paths
                    info = TOML.parsefile(joinpath(path, "package.toml"))
                    msg *= " – $(info["repo"])"
                    break
                end
                msg *= "\n"
            end
            info(msg)
            ambig = true
        end
        ambig && error("interactive package choice not yet implemented")
        Dict(name => first(first(where[name])) for name in names)
    end
    # find applicable package versions
    versions = Dict{String,Dict{VersionNumber,SHA1}}()
    for name in names
        uuid, paths = first(where[name])
        for path in paths
            vers = TOML.parsefile(joinpath(path, "versions.toml"))
            versions[name] = Dict{VersionNumber,SHA1}()
            for (v, d) in vers
                ver = VersionNumber(v)
                ver in pkgs[name] || continue
                versions[name][ver] = SHA1(d["hash-sha1"])
            end
        end
    end
    return versions
end

end # module
