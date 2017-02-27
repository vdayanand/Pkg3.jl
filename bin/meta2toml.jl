#!/usr/bin/env julia

using Base.Random: UUID
using Base.Pkg.Types
using Base.Pkg.Reqs: Reqs, Requirement
using Base: thispatch, thisminor, nextpatch, nextminor
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

struct Require
    versions::VersionInterval
    systems::Vector{Symbol}
end

struct Version
    sha1::String
    julia::VersionInterval
    requires::Dict{String,Require}
end

struct Package
    uuid::UUID
    url::String
    versions::Dict{VersionNumber,Version}
end

function load_requires(path::String)
    requires = Dict{String,Require}()
    isfile(path) || return requires
    for r in filter!(r->r isa Requirement, Reqs.read(path))
        @assert length(r.versions.intervals) == 1
        new = haskey(requires, r.package)
        versions, systems = r.versions.intervals[1], r.system
        if haskey(requires, r.package)
            versions = versions ∩ requires[r.package].versions
            systems  = systems  ∪ requires[r.package].systems
        end
        requires[r.package] = Require(versions, systems)
    end
    return requires
end

function load_versions(dir::String)
    versions = Dict{VersionNumber,Version}()
    isdir(dir) || return versions
    for ver in readdir(dir)
        path = joinpath(dir, ver)
        sha1 = joinpath(path, "sha1")
        isfile(sha1) || continue
        requires = load_requires(joinpath(path, "requires"))
        julia = pop!(requires, "julia", Require(VersionInterval(v"0.1", v"0.6"), []))
        @assert isempty(julia.systems)
        versions[VersionNumber(ver)] = Version(readchomp(sha1), julia.versions, requires)
    end
    return versions
end

function load_packages(dir::String)
    packages = Dict{String,Package}()
    for pkg in readdir(dir)
        path = joinpath(dir, pkg)
        url = joinpath(path, "url")
        versions = joinpath(path, "versions")
        isfile(url) || continue
        packages[pkg] = Package(uuid5(uuid_julia, pkg), readchomp(url), load_versions(versions))
    end
    return packages
end

const julia_versions = [VersionNumber(0,m) for m=1:5]

macro clean(ex) :(x = $(esc(ex)); $(esc(:clean)) &= x; x) end

function prune!(packages::Associative{String,Package})
    while true
        clean = true
        filter!(packages) do p, pkg
            filter!(pkg.versions) do v, ver
                @clean v == thispatch(v) > v"0.0.0" &&
                any(julia in ver.julia for julia in julia_versions) &&
                all(ver.requires) do kv
                    r, req = kv
                    haskey(packages, r) &&
                    any(w->w in req.versions, keys(packages[r].versions))
                end
            end
            @clean !isempty(pkg.versions)
        end
        clean && return packages
    end
end

≲(v::VersionNumber, t::NTuple{0,Int}) = true
≲(v::VersionNumber, t::NTuple{1,Int}) = v.major ≤ t[1]
≲(v::VersionNumber, t::NTuple{2,Int}) = v.major < t[1] ||
                                        v.major ≤ t[1] && v.minor ≤ t[2]
≲(v::VersionNumber, t::NTuple{3,Int}) = v.major < t[1] ||
                                        v.major ≤ t[1] && v.minor < t[2] ||
                                        v.major ≤ t[1] && v.minor ≤ t[2] && v.patch ≤ t[3]

≲(t::NTuple{0,Int}, v::VersionNumber) = true
≲(t::NTuple{1,Int}, v::VersionNumber) = t[1] ≤ v.major
≲(t::NTuple{2,Int}, v::VersionNumber) = t[1] < v.major ||
                                        t[1] ≤ v.major && t[2] ≤ v.minor
≲(t::NTuple{3,Int}, v::VersionNumber) = t[1] < v.major ||
                                        t[1] ≤ v.major && t[2] < v.minor ||
                                        t[1] ≤ v.major && t[2] ≤ v.minor && t[3] ≤ v.patch

function compress_versions(inc::Vector{VersionNumber}, exc::Vector{VersionNumber})
    @assert issorted(inc) && issorted(exc)
    @assert isempty(inc ∩ exc)
    pairs = []
    if isempty(exc)
        a, b = first(inc), last(inc)
        lo = (a.major, a.minor, a.patch)
        hi = (b.major, b.minor, b.patch)
        i, m = 1, min(2, length(lo), length(hi))
        for i = 1:m
            lo[i] == hi[i] || break
        end
        push!(pairs, lo[1:i] => hi[1:i])
    else
        for v in inc
            lo = hi = (v.major, v.minor, v.patch)
            @assert !any(lo ≲ w ≲ hi for w in exc)
            for t in ((v.major, v.minor),)
                !any(t ≲ w ≲ t for w in exc) && (lo = hi = t)
            end
            if isempty(pairs) || any(pairs[end][1] ≲ w ≲ hi for w in exc)
                push!(pairs, lo => hi) # need a new interval
            else
                pairs[end] = pairs[end][1] => hi # can be merged with last
            end
        end
    end
    @assert all(any(p[1] ≲ v ≲ p[2] for p ∈ pairs) for v ∈ inc)
    @assert all(!any(p[1] ≲ v ≲ p[2] for p ∈ pairs) for v ∈ exc)
    return pairs
end

function compress_version_map(fwd::Dict{VersionNumber,X}) where X
    rev = Dict{X,Vector{VersionNumber}}()
    versions = VersionNumber[]
    for (v, x) in fwd
        push!(get!(rev, x, VersionNumber[]), v)
        push!(versions, v)
    end
    sort!(versions)
    Dict(x => compress_versions(sort!(inc), setdiff(versions, inc)) for (x, inc) in rev)
end

version_string(a::Tuple{}, b::Tuple{}) = "*"
version_string(a::NTuple{m,Int}, b::NTuple{n,Int}) where {m,n} =
    a == b ? join(a, '.') : "$(join(a, '.'))-$(join(b, '.'))"
version_string(p::Pair) = version_string(p...)

function print_package_metadata(pkg::String, p::Package)
    print("""
    [$pkg]
    uuid = "$(p.uuid)"
    repository = "$(p.url)"

    """)
end

function print_versions_sha1(pkg::String, p::Package)
    print("""
        [$pkg.versions.sha1]
    """)
    for (ver, v) in sort!(collect(p.versions), by=first)
        print("""
            "$ver" = "$(v.sha1)"
        """)
    end
    println()
end

function compress_julia_versions(julia::VersionInterval)
    inc = filter(v->v in julia, julia_versions)
    exc = setdiff(julia_versions, inc)
    pairs = compress_versions(inc, exc)
    @assert length(pairs) == 1
    return first(pairs)
end

function print_versions_julia(pkg::String, p::Package)
    print("""
        [$pkg.versions.julia]
    """)
    try
    d = Dict(ver => compress_julia_versions(v.julia) for (ver, v) in p.versions)
    for (julias, pairs) in sort!(collect(compress_version_map(d)), by=first∘first∘last)
        lhs = repr(version_string(julias))
        rhs = length(pairs) <= 1 ? repr(version_string(pairs[1])) :
            "[" * join(map(repr∘version_string, pairs), ", ") * "]"
        print("""
            $lhs = $rhs
        """)
    end
    catch
        warn(pkg)
    end
    println()
end

dir = length(ARGS) >= 1 ? ARGS[1] : Pkg.dir("METADATA")
packages = load_packages(dir)
prune!(packages)

if !isinteractive()
    for (pkg, p) in sort!(collect(packages), by=lowercase∘first)
        print_package_metadata(pkg, p)
        print_versions_sha1(pkg, p)
        print_versions_julia(pkg, p)
    end
end
