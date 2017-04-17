## Some graph functions ##

using Iterators

#=
X = sparse([1,1,2,2,3,4,4], [2,5,5,3,4,5,6], 1, 6, 6)
X += X'
G = dropzeros!(1 - X)
=#

#=
G1 = full(sparse(
    [1, 1, 1, 2, 2, 2, 2, 2, 3, 3, 3, 3, 3, 4, 4, 4, 4, 4, 4, 5,
     5, 5, 5, 5, 6, 6, 6, 6, 6, 6, 6, 6, 7, 7, 7, 7, 7, 7, 7, 7,
     8, 8, 8, 8, 8, 9, 9, 9, 9, 9, 10, 10, 10, 10, 11, 11, 11, 11],
    [2, 3, 4, 1, 4, 5, 6, 7, 1, 4, 5, 6, 7, 1, 2, 3, 5, 6, 7, 2,
     3, 4, 6, 7, 2, 3, 4, 5, 8, 9, 10, 11, 2, 3, 4, 5, 8, 9, 10,
     11, 6, 7, 9, 10, 11, 6, 7, 8, 10, 11, 6, 7, 8, 9, 6, 7, 8, 9],
    1
))
p1 = collect(1:Base.LinAlg.checksquare(G1))

G2 = full(sparse(
    [1, 1, 1, 2, 3, 3, 3, 4, 5, 5, 6, 6, 6, 6, 7, 7, 7, 8, 8, 10, 10,
     10, 10, 11, 11, 11, 12, 12, 12, 12, 13, 13, 13, 13, 14, 14, 14],
    [3, 4, 5, 1, 2, 4, 5, 2, 4, 2, 7, 8, 9, 10, 8, 9, 10, 9, 10, 11,
     12, 13, 14, 9, 10, 14, 9, 10, 11, 13, 9, 10, 11, 12, 9, 10, 11],
    1
))
p2 = collect(1:Base.LinAlg.checksquare(G2))
=#

const \ = setdiff

function BronKerboschTomita(emit, G, R, P, X)

    N(G, v) = find(G[:, v])

    let available = unique(pkg_map[v] for V in (R, P) for v in V)
        # if any in R are unsatisfiable, return
        for v in R, r in req_map[v]
            r in available || return true
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
    isempty(P) && isempty(X) && return emit(R) != :break

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
        BronKerboschTomita(emit, G, [R; v], P \ Nv, X \ Nv) || return false
        filter!(x -> x != v, P)
        push!(X, v)
    end
    return true
end

function maximal_indepedent_sets(io::IO, G::AbstractMatrix, inds::Vector{Int}=collect(1:size(G,2)))
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
maximal_indepedent_sets(path::String, G::AbstractMatrix, inds::Vector{Int}=collect(1:size(G,2))) =
    open(io->maximal_indepedent_sets(io, G, inds), path, "w")
maximal_indepedent_sets(G::AbstractMatrix, inds::Vector{Int} = collect(1:size(G,2))) =
    maximal_indepedent_sets(STDOUT, G, inds)

function find_independent_set(G::AbstractMatrix, R::Vector{Int}, inds::Vector{Int}=collect(1:size(G,2)))
    found = Int[]
    G = min.(1, G + I) # make each node its own neighbor
    BronKerboschTomita(G, R, inds \ R, Int[]) do R
        found = sort!(R)
        :break
    end
    return found
end

function satisfiable_pairs(G::AbstractMatrix, inds::Vector{Int}=collect(1:size(G,2)))
    G = min.(1, G + I) # make each node its own neighbor
    S = zeros(Int, length(inds), length(inds))
    for (i, x) in enumerate(inds), (j, y) in enumerate(inds)
        G[x,y] != 0 && continue
        S[i,j] != 0 && continue
        v = find_independent_set(G, [x, y], inds)
        @assert iszero(G[v,v] - I)
        w = findin(inds, v)
        S[w,w] = 1
    end
    return S
end

function is_satisfied(V::Vector{Int})
    provided = unique(pkg_map[v] for v in V)
    required = unique(r for v in V for r in req_map[v])
    required ⊆ provided
end

function is_maximal(G::AbstractMatrix, V::Vector{Int}, inds::Vector{Int} = 1:n)
    minimum(sum(G[V, inds\V], 1)) > 0
end

is_module(G::AbstractMatrix, S::Vector{Int}) = !isempty(S) &&
    all(G[i,k] == G[j,k] && G[k,i] == G[k,j] for i in S for j in S for k in indices(G,2)\S)

all_modules(G) = filter!(S->length(S) > 1 && is_module(G, S), collect(subsets(1:size(G,2))))

overlap(A::Vector, B::Vector) = !isempty(A \ B) && !isempty(A ∩ B) && !isempty(B \ A)

function strong_modules(G)
    modules = all_modules(G)
    return filter(A -> all(B -> !overlap(A, B), modules), modules)
end

function is_modular_permutation(G::AbstractMatrix, p::Vector{Int}; modules=strong_modules(G))
    diffs = map(M->diff(findin(p, M)), modules)
    maximum(maximum, diffs) == 1
end

findin_partition(P, S) = sort!(map(x->findfirst(X->x in X, P), S))

function is_modular_partition(G::AbstractMatrix, P::Vector{Vector{Int}}; modules=strong_modules(G))
    sort!(vcat(P...)) == collect(1:size(G,2)) || error("not a partition")
    diffs = map(M->diff(findin_partition(P, M)), modules)
    maximum(maximum, diffs) == 1
end

## McConnell & Montgolfier 2004: "Linear-time modular decomposition of directed graphs"

function tournament_factorizing_permutation(T::AbstractMatrix)
    x → y = T[y,x] < T[x,y]

    function tfp!(lo::Int, hi::Int)
        lo + 1 < hi || return
        v = p[lo]
        @assert all(v < p[k] for k = lo+1:hi)
        i, j = lo+1, hi
        while true
            while i < j && (p[i] → v); i += 1; end
            while i < j && (v → p[j]); j -= 1; end
            i < j || break
            p[i], p[j] = p[j], p[i]
            @assert (p[i] → v) && (v → p[j])
        end
        @assert i == j
        i -= (v → p[j])
        @assert i <= lo || (p[i] → v)
        @assert i >= hi || (v → p[i+1])
        p[i], p[lo] = v, p[i]
        @assert all(p[k] → v for k = lo:i-1)
        @assert all(v → p[k] for k = i+1:hi)
        tfp!(lo, i-1)
        tfp!(i+1, hi)
    end

    p = collect(1:Base.LinAlg.checksquare(T))
    tfp!(1, length(p))
    return p
end

false &&
for _ = 1:1000
    global n, T, p
    n = rand(3:10)
    T = Int[i != j && rand() < 0.5 for i = 1:n, j = 1:n]
    T .= T .⊻ T' .⊻ tril(ones(Int,n,n),-1)
    @assert T + T' + I == ones(n,n)
    p = tournament_factorizing_permutation(T)
    modules = all_modules(T)
    @assert is_modular_permutation(G, p, modules=modules)
end

## Habib, Paul & Viennot: "Partition refinement techniques: an interesting algorithmic tool kit"

function graph_factorizing_permutation(G::AbstractMatrix, V::Vector{Int}=collect(1:Base.LinAlg.checksquare(G)))

    P = [V]
    center::Int = 0
    pivots::Vector{Vector{Int}} = []
    modules::Vector{Vector{Int}} = []
    first_pivot = Dict{Vector{Int},Int}()

    N_adj(x, X=V) = [y for y in X if y != x && G[x,y] != 0]
    N_non(x, X=V) = [y for y in X if y != x && G[x,y] == 0]

    smaller_larger(A, B) = length(A) <= length(B) ? (A, B) : (B, A)

    function refine!(P, S, x)
        i, between = 0, false
        while (i += 1) <= length(P)
            X = P[i]
            if center in X || x in X
                between = !between
                continue
            end
            Xₐ = X ∩ S
            isempty(Xₐ) && continue
            X = X \ Xₐ
            isempty(X) && continue
            P[i] = X
            insert!(P, i + !between, Xₐ)
            add_pivot(X, Xₐ)
            i += 1
        end
    end

    function add_pivot(X, Xₐ)
        if X in pivots
            push!(pivots, Xₐ)
        else
            S, L = smaller_larger(X, Xₐ)
            push!(pivots, S)
            i = findfirst(modules, X)
            if 0 < i
                modules[i] = L
            else
                push!(modules, L)
            end
        end
    end

    function partition_refinement!(P)
        while init_partition!(P)
            while !isempty(pivots)
                E = pop!(pivots)
                for x in E
                    S = N_adj(x) \ E
                    refine!(P, S, x)
                end
            end
        end
    end

    function init_partition!(P)
        maximum(length, P) <= 1 && return false
        if isempty(modules)
            for (i, X) in enumerate(P)
                length(X) > 1 || continue
                x = get(first_pivot, X, first(X))
                N, A = N_non(x, X), N_adj(x, X)
                splice!(P, i, filter(!isempty, [N, [x], A]))
                S, L = smaller_larger(N, A)
                center = x
                push!(pivots, S)
                push!(modules, L)
                break
            end
        else
            X = shift!(modules)
            x = first(X)
            push!(pivots, [x])
            first_pivot[X] = x
        end
        return true
    end

    partition_refinement!(P)
    return map(first, P)
end

## Capelle, Habib & de Montgolfier 2002: "Graph decompositions and factorizing permutations"

isdefined(:StrongModuleTree) ||
struct StrongModuleTree{T} <: AbstractVector{T}
    kind::Symbol
    edge::Tuple
    nodes::Vector{Union{T,StrongModuleTree{T}}}
end

Base.size(t::StrongModuleTree) = (length(t),)
Base.length(t::StrongModuleTree) = length(t.nodes)
Base.eltype(t::StrongModuleTree) = eltype(t.nodes)
Base.getindex(t::StrongModuleTree, i::Int) = t.nodes[i]

leaf_count(t::StrongModuleTree) = sum(leaf_count, t.nodes)
leaf_count(v::Vector) = sum(leaf_count, v)
leaf_count(x::Any) = 1

first_leaf(t::StrongModuleTree) = first_leaf(first(t.nodes))
first_leaf(v::Vector) = first_leaf(first(v))
first_leaf(x::Any) = x

last_leaf(t::StrongModuleTree) = first_leaf(last(t.nodes))
last_leaf(v::Vector) = last_leaf(last(v))
last_leaf(x::Any) = x

function Base.summary(t::StrongModuleTree)
    edge = t.kind == :prime ? "" : "/$(join(t.edge,"-"))"
    "$(length(t))-node $(t.kind)$edge $(typeof(t))"
end

function Base.show(io::IO, t::StrongModuleTree)
    parens = t.kind == :prime ? "{}" : t.kind == :linear ? "[]" : "()"
    compact = get(io, :compact, false)
    if compact
        print(io, parens[1], first_leaf(t), " × ", leaf_count(t), parens[2])
    else
        print(io, parens[1])
        for (i, x) in enumerate(t)
            print(io, x)
            i < length(t) && print(io, " ")
        end
        print(io, parens[2])
    end
end

Base.getindex(v::Vector{T}, t::StrongModuleTree) where {T} =
    StrongModuleTree{T}(t.kind, t.edge, map(x->v[x], t.nodes))

function StrongModuleTree(G::AbstractMatrix, v::AbstractVector{T}, op::Vector{Int}, cl::Vector{Int}) where T

    function classify_nodes(t::Vector)
        n = length(t)
        counts = zeros(Int, n)
        x, y = first_leaf(t[1]), first_leaf(t[2])
        edge = (G[y,x], G[x,y])
        for i = 1:n, j = 1:n
            i == j && continue
            x, y = first_leaf(t[i]), first_leaf(t[j])
            a, b = G[y,x], G[x,y]
            if edge == (a, b)
                counts[i] += 1
            elseif edge == (b, a)
                counts[j] += 1
            else
                break
            end
        end
        sort!(counts)
        kind = a == b && all(c -> c == n-1, counts) ? :complete :
            all(d -> d == 2, diff(counts)) ? :linear : :prime
        edge[1] <= edge[2] || (edge = reverse(edge))
        kind == :prime && (edge = ())
        StrongModuleTree{T}(kind, edge, map(x->x isa Vector ? classify_nodes(x) : x, t))
    end

    function delete_weak_modules!(t::StrongModuleTree)
        i = 0
        while (i += 1) <= length(t)
            x = t[i]
            x isa StrongModuleTree || continue
            delete_weak_modules!(x)
            t.kind == x.kind && t.edge == x.edge || continue
            splice!(t.nodes, i, x.nodes)
            i += length(x)
        end
    end

    function sort_nodes!(t::StrongModuleTree)
        sort!(t.nodes, by=sort_nodes!)
        return first_leaf(t)
    end
    sort_nodes!(x::Any) = x

    s = Any[[]]
    for (j, x) = enumerate(v)
        for _ = 1:op[j]
            t = []
            push!(s[end], t)
            push!(s, t)
        end
        push!(s[end], x)
        for _ = 1:cl[j]
            pop!(s)
        end
    end
    t = classify_nodes(s[end])
    delete_weak_modules!(t)
    sort_nodes!(t)
    return t
end

function StrongModuleTree(G::AbstractMatrix, p::Vector{Int}=graph_factorizing_permutation(G))
    n = length(p)
    op = zeros(Int,n); op[1] = 1
    cl = zeros(Int,n); cl[n] = 1
    lc = collect(1:n-1)
    uc = collect(2:n)
    # count open and close parens in fracture tree
    # find lower and upper cutters for node pairs
    for j = 1:n-1
        for i = 1:j-1
            G[p[i],p[j]] == G[p[i],p[j+1]] &&
            G[p[j],p[i]] == G[p[j+1],p[i]] && continue
            op[i] += 1
            cl[j] += 1
            lc[j] = i
            break
        end
        j += 1
        for i = n:-1:j+1
            G[p[i],p[j-1]] == G[p[i],p[j]] &&
            G[p[j-1],p[i]] == G[p[j],p[i]] && continue
            op[j] += 1
            cl[i] += 1
            uc[j-1] = i
            break
        end
    end
    # remove non-module "dummy" nodes
    let s = Int[]
        for j = 1:n
            for _ = 1:op[j]; push!(s, j); end
            for _ = 1:cl[j]
                i = pop!(s)
                if i < j
                    l = minimum(lc[k] for k = i:j-1)
                    u = maximum(uc[k] for k = i:j-1)
                    i <= l && u <= j && continue
                end
                op[i] -= 1
                cl[j] -= 1
            end
        end
    end
    # create nodes for consecutive twins
    let s = Int[], t = Int[]
        l = 1
        for k = 1:n
            for _ = 1:op[k]+1
                push!(s, k) # matching node stack
                push!(t, l) # matching twin stack
                l = k
            end
            for c = cl[k]:-1:0
                i = pop!(t)
                j = pop!(s)
                l = i # continue twin chain by default
                i < j || continue
                if i <= lc[j-1] < uc[j-1] <= k
                    # this node and prev are twins
                    if c > 0
                        # not last parens ∴ last twin
                        op[i] += 1
                        cl[k] += 1
                        l = k + 1
                    end
                else # this node and prev not twins
                    if i < j-1
                        op[i] += 1
                        cl[j-1] += 1
                    end
                    l = j # this node starts new chain
                end
            end
        end
    end
    # remove singleton "dummy" nodes
    let s = Int[]
        for j = 1:n
            for _ = 1:op[j]; push!(s, j); end
            i′ = 0
            for _ = 1:cl[j]
                i = pop!(s)
                if i == i′
                    op[i] -= 1
                    cl[j] -= 1
                end
                i′ = i
            end
        end
    end
    op[1] -= 1
    cl[n] -= 1
    # construct and normalize the tree
    return StrongModuleTree(G, p, op, cl)
end

false &&
for _ = 1:1000
    global n, G, p, T, p′, T′
    n = rand(3:10)
    G = Int[i != j && rand() < 0.5 for i = 1:n, j = 1:n]
    G .= G .⊻ G'
    @assert G == G'
    p = graph_factorizing_permutation(G)
    @assert is_modular_permutation(G, p)
    T = StrongModuleTree(G, p)
    for _ = 1:10
        p′ = graph_factorizing_permutation(G, shuffle(1:n))
        @assert is_modular_permutation(G, p′)
        T′ = StrongModuleTree(G, p′)
        @assert T == T′
    end
end
