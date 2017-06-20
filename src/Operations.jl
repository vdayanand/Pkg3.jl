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

const line_re = r"""
    ^ \s*
    ([0-9a-f]{8}-[0-9a-f]{4}-[0-9a-f]{4}-[0-9a-f]{4}-[0-9a-f]{12})
    \s* = \s* \{
    \s* name \s* = \s* "([^"]*)" \s*,
    \s* path \s* = \s* "([^"]*)" \s*,?
    \s* \} \s* $
"""x

function find_registered(names::Vector{String})
    # name --> uuid --> paths
    where = Dict{String,Dict{UUID,Vector{String}}}()
    isempty(names) && return where
    names_re = let p = "\\bname\\s*=\\s*\"(?:$(names[1])"
        for i in 2:length(names)
            p *= "|$(names[i])"
        end
        p *= ")\""
        Regex(p)
    end
    for registry in registries()
        packages_file = joinpath(registry, "packages.toml")
        open(packages_file) do io
            for line in eachline(io)
                ismatch(names_re, line) || continue
                m = match(line_re, line)
                m == nothing && error("misformated packages.toml file: $line")
                uuid = UUID(m.captures[1])
                name = Base.unescape_string(m.captures[2])
                path = Base.unescape_string(m.captures[3])
                if !haskey(where, name)
                    where[name] = Dict{UUID,Vector{String}}()
                end
                if !haskey(where[name], uuid)
                    where[name][uuid] = String[]
                end
                push!(where[name][uuid], abspath(registry, path))
            end
        end
    end
    return where
end

include("libgit2_discover.jl")

default_env() = get(ENV, "JULIA_ENV", nothing)
find_env() = find_env(default_env())

const config_names = ["JuliaConfig.toml", "Config.toml"]
const default_envs = [
    "v$(VERSION.major).$(VERSION.minor).$(VERSION.patch)",
    "v$(VERSION.major).$(VERSION.minor)",
    "v$(VERSION.major)",
    "default",
]

function find_project_env(start_path::String = pwd())
    path = LibGit2.discover(start_path, ceiling = homedir())
    repo = LibGit2.GitRepo(path)
    work = LibGit2.workdir(repo)
    for name in config_names
        path = abspath(work, name)
        isfile(path) && return path
    end
    return abspath(work, config_names[end])
end

function find_default_env()
    for depot in depots(), env in default_envs, name in config_names
        path = joinpath(depot, "environments", env, name)
        isfile(path) && return path
    end
    env = VERSION.major == 0 ? default_envs[2] : default_envs[3]
    return joinpath(user_depot(), "environments", env, config_names[end])
end

function find_env(env::String)
    if isempty(env)
        error("invalid environment name: \"\"")
    elseif env == "/"
        return find_default_env()
    elseif env == "."
        return find_project_env()
    elseif startswith(env, "/") || startswith(env, "./")
        # path to config file or project directory
        splitext(env)[2] == ".toml" && return abspath(env)
        for name in config_names
            path = abspath(env, name)
            isfile(path) && return path
        end
        return abspath(env, config_names[end])
    else # named environment
        for depot in depots()
            path = joinpath(depot, "environments", env, "Config.toml")
            isfile(path) && return path
        end
        return joinpath(user_depot(), "environments", env, "Config.toml")
    end
end

function find_env(::Void)
    try
        return find_project_env()
    catch err
        err isa LibGit2.GitError && err.code == LibGit2.Error.ENOTFOUND || rethrow(err)
    end
    return find_default_env()
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

function load_manifest(env::Union{Void,String}=default_env())
    manifest = find_manifest(env)
    T = Dict{String,Dict{String,Dict{String,String}}}
    return isfile(manifest) ? convert(T, TOML.parsefile(manifest)) : T()
end

function add(pkgs::Dict{String,<:Union{VersionNumber,VersionSpec}})
    names = sort!(collect(keys(pkgs)))
    where = find_registered(names)
    # check for ambiguous package names
    ambig = false
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
    uuids = Dict(name => first(first(where[name])) for name in names)
    # find applicable package versions
    versions = Dict{String,Dict{VersionNumber,SHA1}}()
    for name in names
        uuid, paths = first(where[name])
        for path in paths
            vers = TOML.parsefile(joinpath(path, "versions.toml"))
            versions[name] = Dict{VersionNumber,SHA1}()
            for (v, d) in vers
                ver, spec = VersionNumber(v), pkgs[name]
                (spec isa VersionNumber ? ver == spec : ver in spec) || continue
                versions[name][ver] = SHA1(d["hash-sha1"])
            end
        end
    end
    return versions
end

end # module
