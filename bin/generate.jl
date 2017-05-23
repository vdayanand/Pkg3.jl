#!/usr/bin/env julia

function write_toml(f::Function, names::String...)
    dir = joinpath(names[1:end-1]...)
    mkpath(dir)
    open(joinpath(dir, "$(names[end]).toml"), "w") do io
        f(io, [names...])
    end
end

toml_key(str::String) = ismatch(r"[^\w-]", str) ? repr(str) : str

function toml_section(io::IO, names::Vector{String})
    print(io, '[')
    join(io, map(toml_key, names), '.')
    println(io, ']')
end

cd("tmp") do
    # TOML.parsefile takes ~100 ms for each of these
    # (could be improved with optimizations to TOML)
    write_toml("packages", "uuids") do io, names
        toml_section(io, names)
        for (pkg, p) in sort!(collect(pkgs), by=lowercase∘first)
            println(io, toml_key(pkg), " = ", repr(string(p.uuid)))
        end
    end
    write_toml("packages", "names") do io, names
        toml_section(io, names)
        for (pkg, p) in sort!(collect(pkgs), by=(p->p.uuid.value)∘last)
            println(io, p.uuid, " = ", repr(pkg))
        end
    end
    # bucket package metadata files by first letter
    # TOML.parsefile takes < 20 ms per each file
    buckets = Dict()
    for (pkg, p) in pkgs
        bucket = string(uppercase(first(pkg)))
        push!(get!(buckets, bucket, []), (pkg, p))
    end
    for (bucket, px) in buckets
        write_toml("packages", "metadata", bucket) do io, names
            for (pkg, p) in sort!(px, by=lowercase∘first)
                names[end] = pkg
                toml_section(io, names)
                println(io, "uuid = ", repr(string(p.uuid)))
                println(io, "repo = ", repr(p.url))
                println(io)
            end
        end
    end
end
