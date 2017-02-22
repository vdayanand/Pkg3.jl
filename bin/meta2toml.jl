#!/usr/bin/env julia

using Base.Random: UUID
using Base.Pkg.Types
using Base.Pkg.Reqs: Reqs, Requirement
using Base: thispatch, thisminor, nextpatch, nextminor
using DataStructures: OrderedDict
using SHA

## Computing UUID5 values from (namespace, key) pairs ##

function uuid5(namespace::UUID, key::String)
    data = [reinterpret(UInt8, [namespace.value]); Vector{UInt8}(key)]
    u = reinterpret(UInt128, sha1(data)[1:16])[1]
    u &= 0xffffffffffff0fff3fffffffffffffff
    u |= 0x00000000000050008000000000000000
    return UUID(u)
end
uuid5(namespace::UUID, key::AbstractString) = uuid5(namespace, String(key))

const uuid_dns = UUID(0x6ba7b810_9dad_11d1_80b4_00c04fd430c8)
const uuid_julia = uuid5(uuid_dns, "julialang.org")

## Loading data into various data structures ##

const julia_range = VersionInterval(v"0.1", v"0.5")

struct Require
    versions::VersionInterval
    systems::Vector{Symbol}
end

struct Version
    sha1::String
    julia::VersionInterval
    requires::OrderedDict{String,Require}
    Version(sha1, julia, requires) = new(sha1, julia_range ∩ julia, requires)
end

struct Package
    uuid::UUID
    url::String
    versions::OrderedDict{VersionNumber,Version}
end

function load_requires(path::String)
    requires = OrderedDict{String,Require}()
    isfile(path) || return requires
    reqs = filter!(r->r isa Requirement, Reqs.read(path))
    for r in sort!(reqs, by=r->r.package)
        @assert length(r.versions.intervals) == 1
        requires[r.package] = !haskey(requires, r.package) ?
            Require(r.versions.intervals[1], r.system) :
            Require(
                requires[r.package].versions ∩ r.versions.intervals[1],
                requires[r.package].systems  ∪ r.system,
            )
    end
    return requires
end

function load_versions(dir::String)
    versions = OrderedDict{VersionNumber,Version}()
    isdir(dir) || return versions
    for ver in sort!(map(VersionNumber, readdir(dir)))
        path = joinpath(dir, string(ver))
        sha1 = joinpath(path, "sha1")
        isfile(sha1) || continue
        requires = load_requires(joinpath(path, "requires"))
        julia = pop!(requires, "julia", Require(VersionInterval(),[]))
        @assert isempty(julia.systems)
        versions[ver] = Version(readchomp(sha1), julia.versions, requires)
    end
    return versions
end

function load_packages(dir::String)
    packages = OrderedDict{String,Package}()
    for pkg in sort!(readdir(dir), by=lowercase)
        path = joinpath(dir, pkg)
        url = joinpath(path, "url")
        versions = joinpath(path, "versions")
        isfile(url) || continue
        packages[pkg] = Package(uuid5(uuid_julia, pkg), readchomp(url), load_versions(versions))
    end
    return packages
end

macro clean(ex) :(x = $(esc(ex)); $(esc(:clean)) &= x; x) end

function prune!(packages::Associative{String,Package})
    while true
        clean = true
        filter!(packages) do p, pkg
            filter!(pkg.versions) do v, ver
                @clean thispatch(v) > v"0.0.0" &&
                all(ver.requires) do kv
                    r, req = kv
                    haskey(packages, r) && any(keys(packages[r].versions)) do v′
                        v′ in req.versions
                    end
                end
            end
            @clean !isempty(pkg.versions)
        end
        clean && return packages
    end
end

dir = length(ARGS) >= 1 ? ARGS[1] : Pkg.dir("METADATA")
packages = load_packages(dir)
prune!(packages)


function version_range(vi::VersionInterval, vs::Vector{VersionNumber})
    bef = filter(v->v < vi.lower,  vs)
    inc = filter(v->v in vi,       vs)
    aft = filter(v->vi.upper <= v, vs)
    isempty(inc) && return ""
    lo, hi = first(inc), last(inc)
    shortlo = isempty(bef) || last(bef) < thisminor(lo)
    shorthi = isempty(aft) || nextminor(hi) <= first(aft)
    if shortlo
        slo = "$(lo.major).$(lo.minor)"
        shorthi || (slo *= ".*")
    else
        @assert last(bef) < thispatch(lo)
        slo = "$(thispatch(lo))"
    end
    if shorthi
        shi = "$(hi.major).$(hi.minor)"
        shortlo || (shi *= ".*")
    else
        @assert nextpatch(hi) <= first(aft)
        shi = "$(thispatch(hi))"
    end
    return slo == shi ? slo : "$slo-$shi"
end

macro preamble()
    quote
        path = joinpath(dir, pkg)
        urlf = joinpath(path, "url")
        vers = joinpath(path, "versions")
        isfile(urlf) && isdir(vers) || continue
        versions = sort!(map(VersionNumber, readdir(vers)))
        thispatch(versions[1]) == v"0.0.0" && shift!(versions)
        isempty(versions) && continue
    end |> esc
end

for (i, pkg) in enumerate(names)
    @preamble
    vrange[pkg] = versions
end

for (i, pkg) in enumerate(names)
    @preamble
    # emit package entry
    println("""
    [$pkg]
    uuid = "$(uuid5(uuid_julia, pkg))"
    repository = "$(readchomp(urlf))"
    """)
    for (j, v) in enumerate(versions)
        # emit version info
        verd = joinpath(vers, string(v))
        sha1f = joinpath(verd, "sha1")
        isfile(sha1f) || continue
        reqf = joinpath(verd, "requires")
        reqs = filter!(r->isa(r, Requirement), Reqs.read(reqf))
        reqd, sysd = Dict(), Dict()
        for r in reqs
            p, vs = r.package, r.versions
            reqd[p] = haskey(reqd,p) ? intersect(reqd[p], r.versions) : r.versions
            append!(get!(sysd, p, String[]), r.system)
        end
        julias = pop!(reqd, "julia", VersionSet())
        @assert length(julias.intervals) == 1
        julia = julias.intervals[1]
        jlo = thisminor(max(julia.lower, v"0.1"))
        jhi = thisminor(min(julia.upper, v"0.5"))
        jver = jlo == jhi ? "$(jlo.major).$(jlo.minor)" :
            "$(jlo.major).$(jlo.minor)-$(jhi.major).$(jhi.minor)"
        all(p->haskey(vrange, p), keys(reqd)) || continue
        println("""
            [[$pkg.version]]
            version = "$v"
            julia = "$jver"
            SHA1 = "$(readchomp(sha1f))"
        """)
        # emit compatibility info
        for p in sort!(collect(keys(reqd)))
            vs = reqd[p]
            @assert length(vs.intervals) == 1
            print("""
                    [$pkg.version.$p]
                    uuid = "$(uuid5(uuid_julia, p))"
                    versions = "$(version_range(vs.intervals[1], vrange[p]))"
            """)
            sysp = sysd[p]
            if length(sysp) > 0
                systems = length(sysp) == 1 ? "\"$(sysp[1])\"" :
                    """["$(join(unique(sort!(sysp)), "\", \""))"]"""
                print("""
                        os = $systems
                """)
            end
            println()
        end
        # break
    end
    # break
end
