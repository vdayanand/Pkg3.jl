#!/usr/bin/env julia

using PicoSAT

function representative_versions(pkg::String)
    p = pkgs[pkg]
    versions = sort!(collect(keys(p.versions)))
    rep = pop!(versions)
    reps = [rep]
    while !isempty(versions)
        rep′ = pop!(versions)
        for (pkg′, p′) in pkgs, (ver, v) in p′.versions
            haskey(v.requires, pkg) || continue
            vers = v.requires[pkg].versions
            if (rep in vers) ⊻ (rep′ in vers)
                push!(reps, rep′)
                rep = rep′
                break
            end
        end
    end
    return reverse!(reps)
end

const packages = sort!(collect(keys(pkgs)), by=lowercase)
const versions = Tuple{String,VersionNumber}[]
for pkg in packages, ver in representative_versions(pkg)
    push!(versions, (pkg, ver))
end

function version_to_package()
    n = length(versions)
    pkgs = zeros(Int, n)
    for (i, p1) in enumerate(packages),
        (j, (p2, _)) in enumerate(versions)
        if p1 == p2
            pkgs[j] = i
        end
    end
    return pkgs
end
const ver_to_pkg = version_to_package()

function package_to_version_range()
    m, n = length(packages), length(versions)
    los = zeros(Int, m)
    his = zeros(Int, m)
    for v in 1:n
        p = ver_to_pkg[v]
        los[p] == 0 && (los[p] = v)
        his[p] = v
    end
    [lo:hi for (lo,hi) in zip(los,his)]
end
const pkg_to_vers = package_to_version_range()

function version_to_required_packages()
    ind = Dict(p => i for (i, p) in enumerate(packages))
    [sort!([ind[dep]
        for dep in keys(pkgs[pkg].versions[ver].requires)])
        for (pkg, ver) in versions]
end
const ver_to_reqs = version_to_required_packages()

function incompatibility_matrix()
    n = length(versions)
    X = spzeros(Int, n, n)
    for (i, (p1, v1)) in enumerate(versions)
        r = pkgs[p1].versions[v1].requires
        for (j, (p2, v2)) in enumerate(versions)
            if haskey(r, p2) && v2 ∉ r[p2].versions
                X[i,j] = X[j,i] = 1
            end
        end
    end
    @assert iszero(diag(X))
    @assert issymmetric(X)
    return X
end
const X = incompatibility_matrix()

function build_cnf!(X::AbstractMatrix, cnf::Vector{Vector{Int}})
    n = length(versions)
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

function unsatisfiable(X::AbstractMatrix)
    unsat = Int[]
    cnf = build_cnf!(X, [[0]])
    for v in 1:length(versions)
        print((v, versions[v]), " ...")
        cnf[1][1] = v
        if PicoSAT.solve(cnf) == :unsatisfiable
            print(" UNSAT")
            push!(unsat, v)
        end
        println()
    end
    return unsat
end

function duplicates()
    dups = Int[]
    for p in 1:length(packages)
        d = Dict(X[:,v] => v for v in pkg_to_vers[p])
        append!(dups, setdiff(pkg_to_vers[p], values(d)))
    end
    return dups
end

deleteat!(versions, sort!(unsatisfiable(X) ∪ duplicates()))

append!(empty!(packages), unique(first.(versions)))
append!(empty!(ver_to_pkg), version_to_package())
append!(empty!(pkg_to_vers), package_to_version_range())
append!(empty!(ver_to_reqs), version_to_required_packages())
const X = incompatibility_matrix()

const versions_rev = Dict(v => i for (i, v) in enumerate(versions))
const packages_rev = Dict(p => i for (i, p) in enumerate(packages))

const m, n = length(packages), length(versions)
const P = sparse(1:n, ver_to_pkg, 1, n, m)

@assert all(x->x == 1, sum(P, 2))
@assert all(x->x >= 1, sum(P, 1))

function requirements_matrix()
    pairs = [(i, j) for (j, v) in enumerate(ver_to_reqs) for i in v]
    R = sparse(first.(pairs), last.(pairs), 1, m, n)
    @assert iszero(min.(P*R, P*P'))
    return R
end
const R = requirements_matrix()

function package_to_requiring_versions()
    rev = [Int[] for _ = 1:m]
    for v in 1:n
        for req in ver_to_reqs[v]
            push!(rev[req], v)
        end
    end
    return rev
end
const pkg_to_reqd =  package_to_requiring_versions()

function iterate_dependencies(X, P, R)
    Z = min.(1, X .+ P*P')
    D = dropzeros!(max.(0, min.(1, P*R) .- Z))
    for i = 1:typemax(Int)
        n = countnz(D)
        D .= dropzeros!(max.(0, min.(1, D .+ D^2) .- Z))
        countnz(D) <= n && break
    end
    @assert iszero(min.(D, Z))
    return D
end
const D = iterate_dependencies(X, P, R)

function deep_requirements!()
    cnf = build_cnf!(X, [[0], [0]])
    for (r, v) in zip(findn(P'D)...)
        ver_to_pkg[v] == r && continue
        r in ver_to_reqs[v] && continue
        cnf[1][1], cnf[2][1] = v, -r
        print((v, r, versions[v], packages[r]), " ...")
        if PicoSAT.solve(cnf) == :unsatisfiable
            push!(ver_to_reqs[v], r)
            print(" REQ")
        end
        println()
    end
    foreach(sort!, ver_to_reqs)
end

if !isfile("tmp/ver_to_reqs.jls")
    deep_requirements!()
    open("tmp/ver_to_reqs.jls", "w") do f
        serialize(f, ver_to_reqs)
    end
else
    const ver_to_reqs = open(deserialize, "tmp/ver_to_reqs.jls")
end
R .= requirements_matrix()
dropzeros!(R)

function propagate_consistent_conflicts!()
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
propagate_consistent_conflicts!()
D .= iterate_dependencies()
dropzeros!(D)

function is_satisfied(vers::Vector{Int})
    provided = unique(ver_to_pkg[v] for v in vers)
    required = unique(r for v in vers for r in ver_to_reqs[v])
    required ⊆ provided
end

function pairwise_satisfiability(X, P, D)
    # if two conditions are met then we already know
    # that a pair must be satisfiable:
    #   1. no overlap between their potential dependent packages
    #   2. their potential dependent versions have no conflicts
    # in particular, this includes versions with no requirements
    # that don't explicitly conflict
    S = full(D'max.(P*P',X)*D) .== 0
    cnf = build_cnf!(X, [[0], [0]])
    n = checksquare(X)
    p = sortperm(vec(sum(D, 1)), rev=true)
    for a = 1:n-1, b = a+1:n
        i, j = p[a], p[b]
        X[i,j] == 0 || continue
        S[i,j] == 0 || continue
        println((a, b, versions[i], versions[j]))
        cnf[1][1], cnf[2][1] = i, j
        sat = PicoSAT.solve(cnf)
        sat == :unsatisfiable && continue
        vers = filter(x->1 <= x <= n, sat::Vector{Int})
        @assert iszero(X[vers, vers])
        @assert is_satisfied(vers)
        while true
            absent = setdiff(1:m, (ver_to_pkg[k] for k in vers))
            satisfied = setdiff(find(iszero, sum(R[absent,:],1)), vers)
            compatible = satisfied[find(iszero, sum(X[vers,satisfied],1))]
            isempty(compatible) && break
            push!(vers, compatible[indmin(sum(S[compatible,vers],2))])
        end
        sort!(vers)
        @assert iszero(X[vers, vers])
        @assert is_satisfied(vers)
        S[vers,vers] = 1
        println("SAT: +", length(vers), ", ", 100*countnz(S)/length(S), "%")
    end
    return S
end

if !isfile("tmp/satisfiable.jls")
    S = pairwise_satisfiability(X, P, D)
    open("tmp/satisfiable.jls", "w") do f
        serialize(f, S)
    end
else
    S = open(deserialize, "tmp/satisfiable.jls")
end
X .= sparse(ind2sub(size(S), find(iszero, S))..., 1)
D .= iterate_dependencies()
dropzeros!(X)
dropzeros!(D)

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
X = sparse(ind2sub(size(S), find(iszero, S))..., 1)
propagate_requires!(ver_to_reqs, X)
R = let p = [(i, j) for (j, v) in enumerate(ver_to_reqs) for i in v]
    sparse(first.(p), last.(p), 1, m, n)
end
D1 = max.(0, P'R .- X)
D = iterate_dependencies(X1, P, R)
Dp = min.(1, D*P')

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

end

# D1 .= max.(0, D1 .- X)
# D .= max.(0, D .- X)

# @assert iszero(D1 .& D1') # no mutual dependencies
# @assert iszero(D1 .& X) # explict conflicts can't satisfy requirements
# const G = D1 + X
# @assert G .& G' == X
# @assert max.(0, G .- G') == D1
