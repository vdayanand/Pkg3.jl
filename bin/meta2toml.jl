#!/usr/bin/env julia

using Base.Random: UUID
using Base.Pkg.Types
using Base.Pkg.Reqs: Reqs, Requirement
using Base: thispatch, thisminor, nextpatch, nextminor
using SHA

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

const dir = length(ARGS) >= 1 ? ARGS[1] : Pkg.dir("METADATA")
const names = sort!(readdir(dir), by=lowercase)
const vrange = Dict()

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
