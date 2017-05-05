#!/usr/bin/env julia

using PicoSAT

function version_to_package()
    pkgs = zeros(Int, n)
    for (i, p1) in enumerate(packages),
        (j, (p2, _)) in enumerate(versions)
        if p1 == p2
            pkgs[j] = i
        end
    end
    return pkgs
end

function version_to_required_packages()
    ind = Dict(p => i for (i, p) in enumerate(packages))
    [sort!([ind[dep]
        for dep in keys(pkgs[pkg].versions[ver].requires)])
        for (pkg, ver) in versions]
end

function package_to_requiring_versions()
    rev = [Int[] for _ = 1:m]
    for v in 1:n
        for req in ver_to_reqs[v]
            push!(rev[req], v)
        end
    end
    return rev
end

function package_to_version_range()
    los = zeros(Int, m)
    his = zeros(Int, m)
    for v in 1:n
        p = ver_to_pkg[v]
        los[p] == 0 && (los[p] = v)
        his[p] = v
    end
    [lo:hi for (lo,hi) in zip(los,his)]
end

const ver_to_pkg = version_to_package()
const ver_to_reqs = version_to_required_packages()
const pkg_to_reqd =  package_to_requiring_versions()
const pkg_to_vers = package_to_version_range()

function incompatibility_matrix()
    X = spzeros(Int, n, n)
    for (i, (p1, v1)) in enumerate(versions)
        r = pkgs[p1].versions[v1].requires
        for (j, (p2, v2)) in enumerate(versions)
            if haskey(r, p2) && v2 ∉ r[p2].versions
                X[i,j] = X[j,i] = 1
            end
        end
    end
    return X
end

function iterate_dependencies(X, P, R; verbose=false)
    X = max.(0, X - I) # make each node not its own neighbor
    D = max.(0, min.(1, P'R + I) .- X)
    for i = 1:typemax(Int)
        n = nnz(D)
        verbose && println("Density $i: $(n/length(D))")
        D .= max.(0, min.(1, D^2) .- X)
        nnz(D) <= n && break
    end
    return D
end

function is_satisfied(vers::Vector{Int})
    provided = unique(ver_to_pkg[v] for v in vers)
    required = unique(r for v in vers for r in ver_to_reqs[v])
    required ⊆ provided
end

function build_cnf!(X::AbstractMatrix, cnf::Vector{Vector{Int}})
    for (pkg, vers) in enumerate(pkg_to_vers)
        push!(cnf, [-(n+pkg); vers])
        for ver in vers
            push!(cnf, [-ver, n+pkg])
        end
    end
    for (ver, reqs) in enumerate(ver_to_reqs), req in reqs
        push!(cnf, [-ver, n+req])
    end
    for (v1, v2) in zip(findn(X)...)
        push!(cnf, [-v1, -v2])
    end
    return cnf
end
build_cnf(X::AbstractMatrix) = build_cnf!(X, Vector{Int}[])

function deep_requirements!(ver_to_reqs, P = Main.P, R = Main.R, X = Main.X)
    D = min.(1, P'R)
    while true
        D′ = min.(1, D^2)
        D′ == D && break
        D = D′
    end
    cnf = build_cnf!(X, [[0], [0]])
    for (r, v) in zip(findn(P*D)...)
        ver_to_pkg[v] == r && continue
        r in ver_to_reqs[v] && continue
        cnf[1][1], cnf[2][1] = v, -r
        println((v, r, versions[v], packages[r]))
        sat = PicoSAT.solve(cnf)
        sat != :unsatisfiable && continue
        println("REQ")
        push!(ver_to_reqs[v], r)
    end
    foreach(sort!, ver_to_reqs)
    return ver_to_reqs
end

function pairwise_satisfiability(X::AbstractMatrix, D::AbstractMatrix=Main.D)
    n = checksquare(X)
    S = zeros(Bool, n, n)
    cnf = build_cnf!(X, [[0], [0]])
    p = sortperm(vec(sum(D, 1)), rev=true)
    for a = 1:n-1, b = a+1:n
        i, j = p[a], p[b]
        X[i,j] == 0 || continue
        S[i,j] == 0 || continue
        println((a, b, versions[i], versions[j]))
        cnf[1][1], cnf[2][1] = i, j
        sat = PicoSAT.solve(cnf)
        sat == :unsatisfiable && continue
        vers = filter(x -> 1 <= x <= n, sat::Vector{Int})
        @assert iszero(X[vers, vers])
        @assert is_satisfied(vers)
        while true
            absent = (1:m)\(ver_to_pkg[k] for k in vers)
            satisfied = find(iszero, sum(R[absent,:],1))\vers
            compsat = satisfied[find(iszero, sum(X[vers,satisfied],1))]
            isempty(compsat) && break
            sums = sum(S[:,vers], 2)
            sort!(compsat, by = k -> sums[k])
            for k in compsat
                iszero(X[vers,k]) && push!(vers, k)
            end
        end
        sort!(vers)
        @assert iszero(X[vers, vers])
        @assert is_satisfied(vers)
        S[vers,vers] = 1
        println("SAT: +", length(vers), ", ", 100*countnz(S)/length(S), "%")
    end
    return S
end

# if every compatible version of a requirement requires
# something then make the top version require it too

function propagate_requires!(ver_to_reqs, X::AbstractMatrix = spzeros(Int,n,n))
    # v: a package version
    # req: a required package of v
    # r: a version of req
    # req_req: a requirement of r
    # req_reqs: count requirements of r
    req_reqs = Vector{Int}(m)
    dirty = collect(1:n)
    while !isempty(dirty)
        vers = sort(dirty)
        println("DIRTY: ", length(vers))
        empty!(dirty)
        for v in vers
            println(versions[v])
            for req in ver_to_reqs[v]
                req_vers = filter(x->X[v,x] == 0, pkg_to_vers[req])
                isempty(req_vers) && continue
                req_reqs .= 0
                for r in req_vers
                    for req_req in ver_to_reqs[r]
                        req_reqs[req_req] += 1
                    end
                end
                l = length(req_vers)
                for req_req in find(c->c >= l, req_reqs)
                    if !(req_req in ver_to_reqs[v])
                        push!(ver_to_reqs[v], req_req)
                        append!(dirty, setdiff(pkg_to_reqd[ver_to_pkg[v]], dirty))
                    end
                end
            end
        end
    end
    foreach(sort!, ver_to_reqs)
end

function propagate_conflicts!(X)
    # v: a package version
    # req: a required package of v
    # r: a version of req
    dirty = collect(1:n)
    while !isempty(dirty)
        vers = sort(dirty)
        println("DIRTY: ", length(vers))
        empty!(dirty)
        for v in vers
            println(versions[v])
            for req in ver_to_reqs[v]
                conflicts = spzeros(Int, n)
                for r in pkg_to_vers[req]
                    conflicts += X[:,r]
                end
                l = length(pkg_to_vers[req])
                x = find(c->c == l, conflicts)
                if any(X[x,v] .== 0)
                    X[x,v] = X[v,x] = 1
                    append!(dirty, setdiff(pkg_to_reqd[ver_to_pkg[v]], dirty))
                end
            end
        end
    end
end

if false
P = sparse(ver_to_pkg, 1:n, 1, m, n)
R = let p = [(i, j) for (j, v) in enumerate(ver_to_reqs) for i in v]
    sparse(first.(p), last.(p), 1, m, n)
end
X1 = incompatibility_matrix()
propagate_requires!(ver_to_reqs, X1)
propagate_conflicts!(X1)
@assert issymmetric(X1)

D = iterate_dependencies(X1, P, R)
Dp = min.(1, D*P')

# S = pairwise_satisfiability(X)
S = open(deserialize, "tmp/S.jls")
X = sparse(ind_to_sub(size(S), find(iszero, S))..., 1)
propagate_requires!(ver_to_reqs, X)
R = let p = [(i, j) for (j, v) in enumerate(ver_to_reqs) for i in v]
    sparse(first.(p), last.(p), 1, m, n)
end
D1 = max.(0, P'R .- X)
D = iterate_dependencies(X1, P, R)
Dp = min.(1, D*P')
end

include("graphutils.jl")

# G :: 2^(m+n, m+n)
# p, q ∈ 1:m     (packages)
# v, w ∈ m+(1:n) (versions)
# G[p, q] = 0 (no edges between packages)
# G[p, v] = 1 ⟺ v is a version of p
# G[v, p] = 1 ⟺ v requires p
# G[v, w] = 1 ⟺ v and w are incompatible

function sorttree!(T::StrongModuleTree{Int})
    sort!(T, by=first_leaf)
    sort!(T, by=leaf_count)
    sort!(T, by=node_count)
    return T
end

if false
    TX = sorttree!(StrongModuleTree(X))
    VX = versions[TX]

    TD1 = sorttree!(StrongModuleTree(D1))
    VD1 = versions[TD1]

    TD = sorttree!(StrongModuleTree(D))
    VD = versions[TD]

    G = [spzeros(Int, m, m) P; R' X]
    TG = sorttree!(StrongModuleTree(G))
    VG = [packages; versions][TG]

    H = [spzeros(Int, m, m) P; spzeros(Int, n, m) X]
    TH = sorttree!(StrongModuleTree(H))
    VH = [packages; versions][TH]

    Y = X + D1
    TY = sorttree!(StrongModuleTree(Y))
    VY = versions[TY]
end

if false
    x = find(Dp[:,packages_rev["Colors"]])
    y = find(sum(P[:,x],2))
    Px = P[y,x]
    Rx = R[y,x]
    X1x = X1[x,x]
    pyvx = [packages[y]; versions[x]]

    Xx = ones(Int, length(x), length(x))
    for mis in maximal_indepedent_sets(X1, x)
        z = findin(mis, x)
        Xx[z,z] = 0
    end
    D1x = max.(0, Px'Rx .- Xx)

    Gx = [spzeros(Int, length(y), length(y)) Px; Rx' Xx]
    TG = sorttree!(StrongModuleTree(Gx))
    VG = pyvx[TG]
end

# D1 .= max.(0, D1 .- X)
# D .= max.(0, D .- X)

# @assert iszero(D1 .& D1') # no mutual dependencies
# @assert iszero(D1 .& X) # explict conflicts can't satisfy requirements
# const G = D1 + X
# @assert G .& G' == X
# @assert max.(0, G .- G') == D1