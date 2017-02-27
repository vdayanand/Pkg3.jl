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

macro clean(ex) :(x = $(esc(ex)); $(esc(:clean)) &= x; x) end

function prune!(packages::Associative{String,Package})
    while true
        clean = true
        filter!(packages) do p, pkg
            filter!(pkg.versions) do v, ver
                @clean v == thispatch(v) &&
                thispatch(v) > v"0.0.0" &&
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

version_string(p::Pair) = version_string(p...)
version_string(a::Tuple{}, b::Tuple{}) = "*"

function version_string(a::NTuple{m,Int}, b::NTuple{n,Int}) where {m,n}
    lo, hi = join(a, "."), join(b, ".")
    if a == b
        return lo
    end
    if m + 1 == n && a == b[1:m]
        return b[end] == 0 ? hi : "$hi-"
    end
    if m == n + 1 && a[1:n] == b
        return "$lo+"
    end
    if 0 == m < n
        return all(iszero, b) ? hi : "≤$hi"
    end
    if m > n == 0
        return "≥$lo"
    end
    return "$lo-$hi"
end

function compress_versions(inc::Vector{VersionNumber}, exc::Vector{VersionNumber})
    @assert issorted(inc) && issorted(exc)
    @assert isempty(inc ∩ exc)
    tuples = Tuple[]
    for v in inc
        lo = hi = (v.major, v.minor, v.patch)
        for t in ((v.major, v.minor), (v.major,), ())
            !any(t ≲ w ≤ v for w in exc) && (lo = t)
            !any(v ≤ w ≲ t for w in exc) && (hi = t)
        end
        if isempty(tuples) || any(tuples[end-1] ≲ w ≲ hi for w in exc)
            push!(tuples, lo, hi) # need a new interval
        else
            tuples[end] = hi # can be merged with last
        end
    end
    pairs = [tuples[i] => tuples[i+1] for i = 1:2:length(tuples)]
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
    compressed = Dict{Pair,X}()
    for (x, inc) in rev
        sort!(inc)
        exc = setdiff(versions, inc)
        for r in compress_versions(inc, exc)
            compressed[r] = x
        end
    end
    return compressed
end

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
            $ver = "$(v.sha1)"
        """)
    end
    println()
end

const julia_versions = [VersionNumber(0,m) for m=1:5]

function compress_julia_versions(julia::VersionInterval)
    inc = filter(v->v in julia, julia_versions)
    exc = setdiff(julia_versions, inc)
    compress_versions(inc, exc)
end

function print_versions_julia(pkg::String, p::Package)
    print("""
        [$pkg.versions.julia]
    """)
    d = Dict(ver => compress_julia_versions(v.julia) for (ver, v) in p.versions)
    for (tp, v) in sort!(collect(compress_version_map(d)), by=first∘first)
        lhs = version_string(tp)
        rhs = "\"" * join(map(version_string, v), "\", \"") * "\""
        print("""
            "$lhs" = $rhs
        """)
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
