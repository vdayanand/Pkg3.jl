module Operations

using TOML
using Base.Random: UUID
using Base: LibGit2
using Base: Pkg
using Pkg3.Types

function parse_toml(path::String...; fakeit::Bool=false)
    p = joinpath(path...)
    !fakeit || isfile(p) ? TOML.parsefile(p) : Dict{Any,Any}()
end

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
find_registered(name::String) = find_registered([name])[name]

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
    isfile(config) && return abspath(parse_toml(config)["manifest"])
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
    T = Dict{String,Dict{String,String}}
    return isfile(manifest) ? convert(T, parse_toml(manifest)) : T()
end

inspec(v::VersionNumber, s::VersionSpec) = v in s
inspec(v::VersionNumber, s::VersionNumber) = v == s

function parse_version_set(s::String)
    parts = split(s, '-')
    length(parts) == 1 && return VersionSet(VersionSpec(parts[1]))
    length(parts) != 2 && error("invalid version spec: ", repr(s))
    lower = VersionSet(VersionSpec(parts[1])).intervals[1].lower
    upper = VersionSet(VersionSpec(parts[2])).intervals[1].upper
    return VersionSet(lower, upper)
end

function add(pkgs::Dict{String,<:Union{VersionNumber,VersionSpec}})
    names = sort!(collect(keys(pkgs)))
    regs = find_registered(names)
    manifest = load_manifest()
    uuids = Dict{String,UUID}(name => info["uuid"] for (name, info) in manifest)
    # disambiguate package names
    let ambig = false
        for name in names
            if haskey(uuids, name)
                uuid = uuids[name]
                if haskey(regs, name)
                    for uuid′ in collect(keys(regs[name]))
                        uuid′ == uuid || delete!(regs[name], uuid)
                    end
                end
                haskey(regs, name) && !isempty(regs[name]) && continue
                error("""
                $name/$uuid found in manifest but not in registries;
                To install a different $name package `pkg rm $name` and then do `pkg add $name` again.
                """)
            elseif length(regs[name]) == 1
                uuids[name] = first(first(regs[name]))
            else
                msg = "$name is ambiguous, it could refer to:\n"
                for (i, (uuid, paths)) in enumerate(sort!(collect(regs[name]), by=first))
                    msg *= " [$i] $uuid"
                    for path in paths
                        info = parse_toml(path, "package.toml")
                        msg *= " – $(info["repo"])"
                        break
                    end
                    msg *= "\n"
                end
                info(msg)
                ambig = true
            end
        end
        ambig && error("interactive package choice not yet implemented")
    end
    merge!(regs, find_registered(setdiff(keys(uuids), names)))
    where = Dict(name => paths for (name, info) in regs for (uuid, paths) in info)

    # compute reqs & deps for Pkg.Resolve.resolve
    #  - reqs: what we need to choose versions for
    #  - deps: relevant portion of dependency graph

    for (name, uuid) in uuids
        haskey(pkgs, name) && haskey(manifest, name) || continue
        # if already satisfied ignore add request, otherwise install
        # a requested version and ignore the version in the manifest
        if haskey(manifest[name], "version")
            if inspec(VersionNumber(manifest[name]["version"]), pkgs[name])
                delete!(pkgs, name)
                continue
            end
        end
        delete!(manifest, name)
    end
    # keys of manifest and pkgs are now disjoint

    # reqs :: String --> VersionSet
    reqs = convert(Pkg.Types.Requires, pkgs)

    # deps needs to contain all versions and all potential requirements
    # except those that are already installed at a given version, which
    # can simply be used to filter out the potential candidate versions

    # deps :: String --> VersionNumber --> (SHA1, String --> VersionSet)
    deps = Dict{String,Dict{VersionNumber,Pkg.Types.Available}}()
    names = sort!(collect(keys(reqs)))
    for name in names
        spec = get(pkgs, name, VersionSpec())
        uuid, paths = uuids[name], where[name]
        deps[name] = Dict{VersionNumber,Pkg.Types.Available}()
        for path in paths
            versions = parse_toml(path, "versions.toml")
            requires = parse_toml(path, "requirements.toml", fakeit=true)
            for (v, d) in versions
                ver = VersionNumber(v)
                inspec(ver, spec) || continue
                r = haskey(requires, v) ?
                    Dict(dep => parse_version_set(r) for (dep, r) in requires[v]) :
                    Pkg.Types.Requires()
                if haskey(r, "julia")
                    VERSION in r["julia"] || continue
                    delete!(r, "julia")
                end
                x = deps[name][ver] = Pkg.Types.Available(SHA1(d["hash-sha1"]), r)
                for dep in keys(x.requires)
                    (dep in names || dep == "julia") && continue
                    found = find_registered(dep)
                    @assert length(found) == 1 # TODO: use UUIDs to disambiguate
                    uuids[dep] = first(found)[1]
                    where[dep] = first(found)[2]
                    push!(names, dep)
                end
            end
        end
    end
    deps = Pkg.Query.prune_dependencies(reqs, deps)

    # find applicable package versions
    versions = Dict{UUID,Dict{VersionNumber,SHA1}}()
    for name in names
        uuid, paths = first(where[name])
        for path in paths
            vers = parse_toml(path, "versions.toml")
            versions[uuid] = Dict{VersionNumber,SHA1}()
            for (v, d) in vers
                ver, spec = VersionNumber(v), pkgs[name]
                (spec isa VersionNumber ? ver == spec : ver in spec) || continue
                versions[uuid][ver] = SHA1(d["hash-sha1"])
            end
        end
    end

    return where, versions
end

end # module
