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

prefix = joinpath("registry")

write_toml(prefix, "packages") do io, names
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
        write_toml(prefix, bucket, pkg, "package") do io, names
            println(io, "name = ", repr(pkg))
            println(io, "uuid = ", repr(string(p.uuid)))
            println(io, "repo = ", repr(p.url))
        end
        write_toml(prefix, bucket, pkg, "versions") do io, names
            for (i, (ver, v)) in enumerate(sort!(collect(p.versions), by=first))
                i > 1 && println(io)
                println(io, "[", toml_key(string(ver)), "]")
                println(io, "hash-sha1 = ", repr(v.sha1))
            end
        end
    end
end
