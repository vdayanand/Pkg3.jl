#!/usr/bin/env julia

toml_key(str::String) = ismatch(r"[^\w-]", str) ? repr(str) : str

function write_toml(f::Function, names::String...)
    dir = joinpath(names[1:end-1]...)
    mkpath(dir)
    open(joinpath(dir, "$(names[end]).toml"), "w") do io
        section = join(map(toml_key, names), '.')
        println(io, "[$section]")
        f(io)
    end
end

cd("tmp") do
    write_toml("packages", "uuids") do io
        for (pkg, p) in sort!(collect(pkgs), by=lowercase∘first)
            println(io, "$(toml_key(pkg)) = \"$(p.uuid)\"")
        end
    end
    write_toml("packages", "names") do io
        for (pkg, p) in sort!(collect(pkgs), by=(p->p.uuid.value)∘last)
            println(io, "$(p.uuid) = \"$pkg\"")
        end
    end
    # TOML.parsefile takes 100 ms for each of these
    # (could be improved with optimizations to TOML)
end
