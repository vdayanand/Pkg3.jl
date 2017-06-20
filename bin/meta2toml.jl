#!/usr/bin/env julia

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

## Preprocessing routines ##

function packagelt(a::String, b::String)
    a == "julia" && b != "julia" && return true
    a != "julia" && b == "julia" && return false
    return lowercase(a) < lowercase(b)
end

function compat_julia(p::Package)
    fwd = Dict(ver => compress_versions(v.julia, julia_versions()) for (ver, v) in p.versions)
    rev = Dict(jul => compress_versions(vers, keys(fwd)) for (jul, vers) in invert_map(fwd))
    return sort!(collect(flatten_keys(invert_map(rev))), by=first∘first)
end

function compat_versions(p::Package, packages=Main.packages)
    fwd = Dict{String,Dict{VersionNumber,Any}}()
    for (ver, v) in p.versions, (req, r) in v.requires
        d = get!(fwd, req, Dict{VersionNumber,Any}())
        d[ver] = compress_versions(r.versions, keys(packages[req].versions))
    end
    vers = sort!(collect(keys(p.versions)))
    uniform = Dict{String,Vector{Any}}()
    nonunif = Dict{String,Vector{Pair}}()
    for (req, d) in sort!(collect(fwd), by=lowercase∘first)
        r = Dict(rv => compress_versions(pv, vers) for (rv, pv) in invert_map(d))
        if length(r) == 1 && length(d) == length(vers)
            # same requirements for all pkg versions
            uniform[req] = first(keys(r))
        else
            # different requirements for various versions
            nonunif[req] = sort!(collect(flatten_keys(invert_map(r))), by=first∘first)
        end
    end
    (sort!(collect(uniform), by=first, lt=packagelt),
     sort!(collect(nonunif), by=first, lt=packagelt))
end

## Package info output routines ##

function print_package_metadata(pkg::String, p::Package; julia=compat_julia(p))
    print("""
    [$pkg]
    uuid = "$(p.uuid)"
    repo = "$(p.url)"
    """)
    if length(julia) == 1
        print("""
        julia = $(versions_repr(julia[1][2]))
        """)
    end
    println()
end

function print_versions_sha1(pkg::String, p::Package)
    print("""
    \t[$pkg.versions.sha1]
    """)
    for (ver, v) in sort!(collect(p.versions), by=first)
        print("""
        \t"$ver" = "$(v.sha1)"
        """)
    end
    println()
end

function print_compat_julia(pkg::String, p::Package; julia=compat_julia(p))
    length(julia) == 1 && return
    print("""
    \t[$pkg.julia]
    """)
    for (versions, julias) in julia
        @assert length(julias) == 1
        print("""
        \t$(versions_repr(versions)) = $(versions_repr(julias))
        """)
    end
    println()
end

# NOTE: UUID mapping is optional. When a dependency is used by name with no
# corresponding UUID mapping, then the current meaning of this name in the same
# registry is assumed. If no UUID mappings are present, the section may be
# skipped entirely. If a package has a dependency which is not registered in the
# same registry, then it  must include a UUID mapping entry for that dependency.

# This code doesn't emit this, but if a dependency name is used to refer to
# different packages across different versions of a package, then this is
# expressed by associating the name with a `version => uuid` table instead of
# just a single UUID value.

function print_compat_uuids(pkg::String, p::Package; pkgs=Main.pkgs)
    print("""
    \t[$pkg.compat.uuids]
    """)
    pkgs = Set{String}()
    for (ver, v) in p.versions
        union!(pkgs, keys(v.requires))
    end
    for pkg in sort!(collect(pkgs), by=lowercase)
        print("""
        \t$pkg = "$(pkgs[pkg].uuid)"
        """)
    end
    println()
end

function print_compat_versions(pkg::String, p::Package; pkgs=Main.pkgs)
    uniform, nonunif = compat_versions(p, pkgs)
    if !isempty(uniform)
        print("""
        \t[$pkg.compat.versions]
        """)
        for (req, v) in uniform
            print("""
            \t$req = $(versions_repr(v))
            """)
        end
        println()
    end
    for (req, r) in nonunif
        print("""
        \t[$pkg.compat.versions.$req]
        """)
        for (pv, rv) in sort!(collect(r), by=first∘first)
            print("""
            \t$(versions_repr(pv)) = $(versions_repr(rv))
            """)
        end
        println()
    end
end

function print_compat(pkg::String, p::Package; pkgs=Main.pkgs)
    julia = compat_julia(p)
    uniform, nonunif = compat_versions(p, pkgs)

    if length(julia) == 1 || !isempty(uniform)
        println("[$pkg]")
        if length(julia) == 1
            println("julia = $(versions_repr(julia[1][2]))")
        end
        for (req, v) in uniform
            println("$req = $(versions_repr(v))")
        end
        println()
    end
    if length(julia) != 1
        println("[$pkg.julia]")
        for (versions, julias) in julia
            println("$(versions_repr(versions)) = $(versions_repr(julias))")
        end
        println()
    end
    for (req, r) in nonunif
        println("[$pkg.$req]")
        for (pv, rv) in r
            println("$(versions_repr(pv)) = $(versions_repr(rv))")
        end
        println()
    end
end

if !isinteractive()
    for (pkg, p) in sort!(collect(pkgs), by=lowercase∘first)
        julia = compat_julia(p)
        print_package_metadata(pkg, p, julia=julia)
        print_versions_sha1(pkg, p)
        print_compat_julia(pkg, p, julia=julia)
        # NOTE: because of optional UUID mapping, this section is totally
        # unnecessary while translating metadata. We could however, represent
        # the Stats => StatsBase rename correctly.
        false && print_compat_uuids(pkg, p)
        print_compat_versions(pkg, p)
    end
end
