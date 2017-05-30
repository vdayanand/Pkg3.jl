#!/usr/bin/env julia

function write_toml(f::Function, names::String...)
    path = joinpath(names...) * ".toml"
    mkpath(dirname(path))
    open(path, "w") do io
        f(io, [names...])
    end
end

toml_key(str::String) = ismatch(r"[^\w-]", str) ? repr(str) : str
toml_key(strs::String...) = join(map(toml_key, [strs...]), '.')

prefix = joinpath("tmp", "packages")

write_toml(prefix, "index") do io, names
    for (pkg, p) in sort!(collect(pkgs), by=(p->p.uuid.value)∘last)
        println(io, p.uuid, " = ", repr(pkg))
    end
end

buckets = Dict()
for (pkg, p) in pkgs
    bucket = string(uppercase(first(pkg)))
    push!(get!(buckets, bucket, []), (pkg, p))
end

for (bucket, b_pkgs) in buckets
    for (pkg, p) in b_pkgs
        write_toml(prefix, "metadata", bucket, pkg) do io, names
            println(io, "name = ", repr(pkg))
            println(io, "uuid = ", repr(string(p.uuid)))
            println(io, "repo = ", repr(p.url))
            for (ver, v) in sort!(collect(p.versions), by=first)
                println(io)
                println(io, "[", toml_key("version", string(ver)), "]")
                println(io, "hash-sha1 = ", repr(v.sha1))
            end
        end
    end
end
