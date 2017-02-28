#!/usr/bin/env julia

using Base.Random: UUID
using Base.Pkg.Types
using Base.Pkg.Reqs: Reqs, Requirement
using Base: thispatch, thisminor, nextpatch, nextminor
using SHA

## General utility functions ##

function invert_map(fwd::Dict{K,V}) where {K,V}
    rev = Dict{V,Vector{K}}()
    for (k, v) in fwd
        push!(get!(rev, v, K[]), k)
    end
    foreach(sort!, values(rev))
    return rev
end

function invert_map(fwd::Dict{Vector{K},V}) where {K,V}
    rev = Dict{V,Vector{K}}()
    for (k, v) in fwd
        append!(get!(rev, v, K[]), k)
    end
    foreach(sort!, values(rev))
    return rev
end

flatten_keys(d::Dict{Vector{K},V}) where {K,V} =
    Dict{K,V}(k => v for (ks, v) in d for k in ks)

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

julia_versions(f::Function) = collect(Iterators.filter(f, VersionNumber(0,m) for m=1:5))
julia_versions(vi::VersionInterval) = julia_versions(v->v in vi)
@eval julia_versions() = $(julia_versions(v->true))

macro clean(ex) :(x = $(esc(ex)); $(esc(:clean)) &= x; x) end

function prune!(packages::Associative{String,Package})
    while true
        clean = true
        filter!(packages) do pkg, p
            filter!(p.versions) do ver, v
                @clean ver == thispatch(ver) > v"0.0.0" &&
                !isempty(julia_versions(v.julia)) &&
                all(v.requires) do kv
                    req, r = kv
                    haskey(packages, req) &&
                    any(w->w in r.versions, keys(packages[req].versions))
                end
            end
            @clean !isempty(p.versions)
        end
        clean && return packages
    end
end

## Functions for representing package version info ##

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

function compress_versions(inc::Vector{VersionNumber}, from::Vector{VersionNumber})
    issorted(inc) || (inc = sort(inc))
    exc = sort!(setdiff(from, inc))
    pairs = []
    if isempty(exc)
        lo, hi = first(inc), last(inc)
        push!(pairs, (lo.major, lo.minor) => (hi.major, hi.minor))
    else
        for v in inc
            t = (v.major, v.minor)
            if any(t ≲ w ≲ t for w in exc)
                t = (v.major, v.minor, v.patch)
            end
            if isempty(pairs) || any(pairs[end][1] ≲ w ≲ t for w in exc)
                push!(pairs, t => t) # need a new interval
            else
                pairs[end] = pairs[end][1] => t # can be merged with last
            end
        end
    end
    @assert all(any(p[1] ≲ v ≲ p[2] for p ∈ pairs) for v ∈ inc)
    @assert all(!any(p[1] ≲ v ≲ p[2] for p ∈ pairs) for v ∈ exc)
    return pairs
end
compress_versions(f::Function, from::Vector{VersionNumber}) =
    compress_versions(filter(f, from), from)
compress_versions(vi::VersionInterval, from::Vector{VersionNumber}) =
    compress_versions(v->v in vi, from)
compress_versions(inc, from) = compress_versions(inc, collect(from))

versions_string(p::Pair) = versions_string(p...)
versions_string(a::Tuple{}, b::Tuple{}) = "*"
versions_string(a::NTuple{m,Int}, b::NTuple{n,Int}) where {m,n} =
    a == b ? join(a, '.') : "$(join(a, '.'))-$(join(b, '.'))"

versions_repr(x) = repr(versions_string(x))
versions_repr(v::Vector) = length(v) == 1 ? repr(versions_string(v[1])) :
    "[" * join(map(repr∘versions_string, v), ", ") * "]"

## Package info output routines ##

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

function print_versions_julia(pkg::String, p::Package)
    print("""
        [$pkg.versions.julia]
    """)
    fwd = Dict(ver => compress_versions(v.julia, julia_versions()) for (ver, v) in p.versions)
    rev = Dict(jul => compress_versions(vers, keys(fwd)) for (jul, vers) in invert_map(fwd))
    for (vers, jul) in sort!(collect(flatten_keys(invert_map(rev))), by=first∘first)
        @assert length(jul) == 1
        print("""
            $(versions_repr(vers)) = $(versions_repr(jul))
        """)
    end
    println()
end

## Load package data and generate registry ##

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
