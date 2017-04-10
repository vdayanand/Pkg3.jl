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
    requires::Dict{String,Require}
end

struct Package
    uuid::UUID
    url::String
    versions::Dict{VersionNumber,Version}
end

Require(versions::VersionInterval) = Require(versions, Symbol[])
Version(sha1::String) = Version(sha1, Dict{String,Require}())

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
        versions[VersionNumber(ver)] = Version(readchomp(sha1), requires)
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
        packages[pkg] = Package(
            uuid5(uuid_julia, pkg),
            readchomp(url),
            load_versions(versions),
        )
    end
    packages["julia"] = Package(
        uuid5(uuid_julia, "julia"),
        "https://github.com/JuliaLang/julia.git",
        Dict(
            v"0.1.1"  => Version("bddb70374fdbb2256782398c13430fa4d8b55d0b"),
            v"0.1.2"  => Version("1187040d027210cc466f0b6a6d54118fd692cf2d"),
            v"0.2.0"  => Version("68e4f269049b3a77eee1e185e619c3513034b3a6"),
            v"0.2.1"  => Version("7c5dea64f55afb7d9a7224a337d3c0bd5a707c50"),
            v"0.3.0"  => Version("768187890c2709bf2ff06818f40e1cdc79bd44b0"),
            v"0.3.1"  => Version("c03f413bbdb46c00033f4eaad402995cfe3b7be5"),
            v"0.3.10" => Version("c8ceeefcc1dc25953a644622a895a3adcbc80dad"),
            v"0.3.11" => Version("483dbf5279efe29f94cf35d42de0ab0a673fb34e"),
            v"0.3.12" => Version("80aa77986ec84067f4696439f10fc77c93703bc7"),
            v"0.3.2"  => Version("8227746b95146c2921f83d2ae5f37ecd146592d8"),
            v"0.3.3"  => Version("b24213b893ace6ca601722f6cdaac15f2b034b7b"),
            v"0.3.4"  => Version("33920266908e950936ef3f8503810a2d47333741"),
            v"0.3.5"  => Version("a05f87b79ad62beb033817fdfdefa270c9557aaf"),
            v"0.3.6"  => Version("0c24dca65c031820b91721139f0291068086955c"),
            v"0.3.7"  => Version("cb9bcae93a32b42cec02585c387396ff11836aed"),
            v"0.3.8"  => Version("79599ada444fcc63cb08aed64b4ff6a415bb4d29"),
            v"0.3.9"  => Version("31efe690bea1fd8f4e44692e205fb72d34f50ad3"),
            v"0.4.0"  => Version("0ff703b40afddf9b705bd6a06d3a59cb4c089ea5"),
            v"0.4.1"  => Version("cbe1bee3a8f27ba4f349556857bf4615ee5fa68f"),
            v"0.4.2"  => Version("bb73f3489d837e3339fce2c1aab283d3b2e97a4c"),
            v"0.4.3"  => Version("a2f713dea5ac6320d8dcf2835ac4a37ea751af05"),
            v"0.4.4"  => Version("ae683afcc4cc824dcf68251569a6cf8760e4fb6f"),
            v"0.4.5"  => Version("2ac304dfba75fad148d4070ef4f8a2e400c305bb"),
            v"0.4.6"  => Version("2e358ce975029ec97aba5994c17d4a2169c3b085"),
            v"0.4.7"  => Version("ae26b25d43317d7dd3ca05f60b70677aab9c0e08"),
            v"0.5.0"  => Version("3c9d75391c72d7c32eea75ff187ce77b2d5effc8"),
            v"0.5.1"  => Version("6445c82d0060dbe82b88436f0f4371a4ee64d918"),
        ),
    )
    return packages
end

@eval julia_versions() = $([VersionNumber(0,m) for m=1:5])
julia_versions(f::Function) = filter(f, julia_versions())
julia_versions(vi::VersionInterval) = julia_versions(v->v in vi)

macro clean(ex) :(x = $(esc(ex)); $(esc(:clean)) &= x; x) end

function prune!(packages::Associative{String,Package})
    while true
        clean = true
        filter!(packages) do pkg, p
            filter!(p.versions) do ver, v
                @clean ver == thispatch(ver) > v"0.0.0" &&
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

## Preprocessing routines ##

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
    (sort!(collect(uniform), by=lowercase∘first),
     sort!(collect(nonunif), by=lowercase∘first))
end

## Load package data ##

const dir = length(ARGS) >= 1 ? ARGS[1] : Pkg.dir("METADATA")
const pkgs = load_packages(dir)
prune!(pkgs)
const versions = [(pkg, ver) for (pkg, p) in pkgs for (ver, v) in p.versions]
sort!(versions, by=last)
sort!(versions, by=lowercase∘first)
const packages = unique(first(v) for v in versions)
const packages_rev = Dict(p => i for (i, p) in enumerate(packages))
const m, n = length(pkgs), length(versions)

## Some computational utility functions ##

function package_map()
    pkgs = zeros(Int, n)
    for (i, p1) in enumerate(packages),
        (j, (p2, _)) in enumerate(versions)
        if p1 == p2
            pkgs[j] = i
        end
    end
    return pkgs
end

function requires_map()
    ind = Dict(p => i for (i, p) in enumerate(packages))
    [sort!([ind[dep] for dep in keys(pkgs[pkg].versions[ver].requires)])
                     for (pkg, ver) in versions]
end

function incompatibility_matrix()
    G = spzeros(n, n)
    for (i, (p1, v1)) in enumerate(versions),
        (j, (p2, v2)) in enumerate(versions)
        r = pkgs[p1].versions[v1].requires
        if haskey(r, p2) && v2 ∉ r[p2].versions
            G[i,j] = G[j,i] = 1
        end
    end
    return G
end

function iterate_dependencies(G, P, R; verbose=false)
    G = max.(0, G - I) # make each node not its own neighbor
    D = max.(0, min.(1, P'R + I) .- G)
    for i = 1:typemax(Int)
        n = nnz(D)
        verbose && println("Density $i: $(n/length(D))")
        D .= max.(0, min.(1, D^2) .- G)
        nnz(D) <= n && break
    end
    return D
end

const pkg_map = package_map()
const req_map = requires_map()
const G = incompatibility_matrix()
const P = sparse(pkg_map, 1:n, 1.0, m, n)
const R = let p = [(i, j) for (j, v) in enumerate(req_map) for i in v]
    sparse(first.(p), last.(p), 1.0, m, n)
end
const D = iterate_dependencies(G, P, R)
const Dp = min.(1, D*P')

## Some graph functions ##

#=
X = sparse([1,1,2,2,3,4,4], [2,5,5,3,4,5,6], 1, 6, 6)
X += X'
G = dropzeros!(1 - X)
=#

"Neighborhood of a node in the graph G, represented as a matrix."
N(G, v) = find(G[:, v])
const \ = setdiff

function BronKerboschTomita(emit, G, R, P, X)
    let available = unique(pkg_map[v] for V in (R, P) for v in V)
        # if any in R are unsatisfiable, return
        for v in R, r in req_map[v]
            r in available || return
        end
        # scrub unsatisfiable versions from P & X
        for V in (P, X)
            filter!(V) do v
                all(r in available for r in req_map[v])
            end
        end
    end
    @show length(R), length(P), length(X)
    # recursion base case
    isempty(P) && isempty(X) && (emit(R); return)
    # pivot: u in P ∪ X minimizing P ∩ N(G, u)
    u, m = 0, typemax(Int)
    for V in (P, X), v in V
        n = sum(G[P, v])
        n < m && ((u, m) = (v, n))
    end
    @assert u != 0
    # recursion
    for v in P ∩ N(G, u)
        Nv = N(G, v)
        BronKerboschTomita(emit, G, [R; v], P \ Nv, X \ Nv)
        filter!(x -> x != v, P)
        push!(X, v)
    end
end

function maximal_indepedent_sets(io::IO, G::AbstractMatrix, inds::Vector{Int} = collect(1:size(G,2)))
    G = min.(1, G + I) # make each node its own neighbor
    M = Vector{Vector{Int}}()
    BronKerboschTomita(G, Int[], copy(inds), Int[]) do R
        push!(M, sort!(R))
        join(io, R, ',')
        println(io)
        flush(io)
    end
    return sort!(M, lt=lexless)
end
maximal_indepedent_sets(path::String, G::AbstractMatrix, inds::Vector{Int} = collect(1:size(G,2))) =
    open(io->maximal_indepedent_sets(io, G, inds), path, "w")
maximal_indepedent_sets(G::AbstractMatrix, inds::Vector{Int} = collect(1:size(G,2))) =
    maximal_indepedent_sets(STDOUT, G, inds)

function is_satisfied(V::Vector{Int})
    provided = unique(pkg_map[v] for v in V)
    required = unique(r for v in V for r in req_map[v])
    required ⊆ provided
end

function is_maximal(G::AbstractMatrix, V::Vector{Int}, inds::Vector{Int} = 1:n)
    minimum(sum(G[V, inds\V], 1)) > 0
end

is_module(G::AbstractMatrix, S::Vector{Int}) = !isempty(S) &&
    all(G[i,k] == G[j,k] for i in S for j in S for k in indices(G,2)\S)

overlap(A::Vector, B::Vector) = !isempty(A \ B) && !isempty(A ∩ B) && !isempty(B \ A)

## McConnell & Montgolfier 2004: "Linear-time modular decomposition of directed graphs"

#=
n = 5
G = rand(n, n)
T = float(G .< G')
@assert T + T' + I == ones(T)
p = tournament_factorizing_permutation(T)
modules = filter(S->is_module(T, S), collect(subsets(1:n)))
@assert all(M->all(x->x == 1, diff(findin(M, p))), modules)
=#

tournament_factorizing_permutation(T::AbstractMatrix) =
    tfp!(T, collect(1:Base.LinAlg.checksquare(T)))

function tfp!(T::AbstractMatrix, p::Vector{Int}, lo::Int=1, hi::Int=length(p))
    if hi - lo > 1
        v = p[lo]
        i, j = lo + 1, hi
        while true
            while i < j && T[v,p[i]] == 0; i += 1; end;
            while i < j && T[p[j],v] == 0; j -= 1; end;
            i < j || break
            p[i], p[j] = p[j], p[i]
        end
        p[j], p[lo] = v, p[j]
        tfp!(T, p, lo, j-1)
        tfp!(T, p, j+1, hi)
    end
    return p
end

#=
G = full(sparse(
    [1, 1, 1, 2, 2, 2, 2, 2, 3, 3, 3, 3, 3, 4, 4, 4, 4, 4, 4, 5,
     5, 5, 5, 5, 6, 6, 6, 6, 6, 6, 6, 6, 7, 7, 7, 7, 7, 7, 7, 7,
     8, 8, 8, 8, 8, 9, 9, 9, 9, 9, 10, 10, 10, 10, 11, 11, 11, 11],
    [2, 3, 4, 1, 4, 5, 6, 7, 1, 4, 5, 6, 7, 1, 2, 3, 5, 6, 7, 2,
     3, 4, 6, 7, 2, 3, 4, 5, 8, 9, 10, 11, 2, 3, 4, 5, 8, 9, 10,
     11, 6, 7, 9, 10, 11, 6, 7, 8, 10, 11, 6, 7, 8, 9, 6, 7, 8, 9],
    1.0
))
n = Base.LinAlg.checksquare(G)
=#

#=
n = 6
G = Int[i != j && rand() < 0.5 for i = 1:n, j = 1:n]
G .= G .⊻ G'
@assert G == G'
p = graph_factorizing_permutation(G)
modules = filter!(S->is_module(G, S), collect(subsets(1:n)))
@assert all(M->all(x->x == 1, diff(findin(M, p))), modules)
strong = filter(A -> all(B -> !overlap(A, B), modules), modules)
=#

graph_factorizing_permutation(G::AbstractMatrix) =
    gfp!(G, collect(1:Base.LinAlg.checksquare(G)))

function gfp!(G::AbstractMatrix, p::Vector{Int}, lo::Int=1, hi::Int=length(p))
    if hi - lo > 1
        v = p[lo]
        i, j = lo + 1, hi
        while true
            while i < j && G[v,p[i]] == 0; i += 1; end;
            while i < j && G[v,p[j]] != 0; j -= 1; end;
            i < j || break
            p[i], p[j] = p[j], p[i]
        end
        gfp!(G, p, lo+1, j)
        gfp!(G, p, j+1, hi)
    end
    return p
end

## Uno & Yagiura 2000: "Fast algorithms to enumerate all common intervals of two permutations"

#=
Brute force all common intervals:
[(x,y) for x=1:n-1 for y=x+1:n for l=1:n-1 for u=l+1:n if sort!(p1[x:y]) == sort!(p2[l:u])]
=#

function all_common_intervals(emit::Function, p::Vector{Int})
    for x = 1:length(p)-1
        l = u = p[x]
        for y = x+1:length(p)
            v = p[y]
            l, u = min(v, l), max(v, u)
            y - x < u - l && continue
            y - x > u - l && break
            emit(x, y)
        end
    end
end

function all_common_intervals(p::Vector{Int})
    intervals = NTuple{2,Int}[]
    all_common_intervals(p) do x, y
        push!(intervals, (x, y))
    end
    return intervals
end

all_common_intervals(emit::Function, p1::Vector{Int}, p2::Vector{Int}) =
    all_common_intervals(emit, invperm(p2)[p1])
all_common_intervals(p1::Vector{Int}, p2::Vector{Int}) =
    all_common_intervals(invperm(p2)[p1])

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
