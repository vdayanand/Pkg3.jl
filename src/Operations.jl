module Operations

using TOML
using Base.Random: UUID
using Pkg3.Types

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

function add(pkgs::Dict{String,<:VersionSpec})
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
            versions[name] = Dict(VersionNumber(v) => SHA1(d["hash-sha1"]) for (v,d) in vers)
        end
    end
    return versions
end

end # module
