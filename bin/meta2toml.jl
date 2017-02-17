#!/usr/bin/env julia

using SHA
# using TOML # not the registered one, but https://github.com/StefanKarpinski/TOML.jl
using Base.Random: UUID

function uuid5(namespace::UUID, key::String)
    data = [reinterpret(UInt8, [namespace.value]); Vector{UInt8}(key)]
    u = reinterpret(UInt128, sha1(data)[1:16])[1]
    u &= 0xffffffffffff0fff3fffffffffffffff
    u |= 0x00000000000050008000000000000000
    return UUID(u)
end

const uuid_dns = UUID(0x6ba7b810_9dad_11d1_80b4_00c04fd430c8)
const uuid_julia = uuid5(uuid_dns, "julialang.org")

const dir = length(ARGS) >= 1 ? ARGS[1] : Pkg.dir("METADATA")
const names = sort!(readdir(dir), by=lowercase)
const packages = Dict()

for (i, pkg) in enumerate(names)
    path = joinpath(dir, pkg)
    urlf = joinpath(path, "url")
    isfile(urlf) || continue
    url = readchomp(urlf)
    uuid = uuid5(uuid_julia, pkg)
    print("""
    [$pkg]
    uuid = "$uuid"
    repository = "$url"
    """)
    i < length(names) && println()
end
