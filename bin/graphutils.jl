## Some graph functions ##

using Iterators
using Combinatorics
using Base.LinAlg: checksquare

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
p1 = collect(1:checksquare(G1))

G2 = full(sparse(
    [1, 1, 1, 2, 3, 3, 3, 4, 5, 5, 6, 6, 6, 6, 7, 7, 7, 8, 8, 10, 10,
     10, 10, 11, 11, 11, 12, 12, 12, 12, 13, 13, 13, 13, 14, 14, 14],
    [3, 4, 5, 1, 2, 4, 5, 2, 4, 2, 7, 8, 9, 10, 8, 9, 10, 9, 10, 11,
     12, 13, 14, 9, 10, 14, 9, 10, 11, 13, 9, 10, 11, 12, 9, 10, 11],
    1
))
p2 = collect(1:checksquare(G2))

n = checksquare(G2)
p = shuffle(1:n)
G = G2[invperm(p),invperm(p)]
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
        iszero(G[R,R] - I) || return :continue
        found = sort!(R)
        return :break
    end
    return found
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

function all_modules(G::AbstractMatrix)
    n = checksquare(G)
    filter!(S->1 < length(S) && is_module(G, S), collect(subsets(1:n)))
end

overlap(A::AbstractVector, B::AbstractVector) =
    # TODO: efficient overlap checking for sorted vectors
    # return true as soon as one elt of each is found:
    # x ∩ y, x \ y, y \ x
    A !== B && !isempty(A \ B) && !isempty(A ∩ B) && !isempty(B \ A)

function strong_modules(G::AbstractMatrix)
    modules = all_modules(G)
    return filter(A -> all(B -> !overlap(A, B), modules), modules)
end

function is_modular_permutation(G::AbstractMatrix, p::Vector{Int}; modules=strong_modules(G))
    isempty(modules) && return true
    diffs = map(M->diff(findin(p, M)), modules)
    maximum(maximum, diffs) == 1
end

findin_partition(P, S) = sort!(map(x->findfirst(X->x in X, P), S))

function is_modular_partition(G::AbstractMatrix, P::Vector{Vector{Int}}; modules=strong_modules(G))
    sort!(vcat(P...)) == collect(1:size(G,2)) || error("not a partition")
    diffs = map(M->diff(findin_partition(P, M)), modules)
    maximum(maximum, diffs) == 1
end

## Habib, Paul & Viennot: "Partition refinement techniques: an interesting algorithmic tool kit"

function graph_factorizing_permutation(G::AbstractMatrix, V::Vector{Int}=collect(1:checksquare(G)))

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
            insert!(P, i + between, Xₐ)
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
                A, N = N_adj(x, X), N_non(x, X)
                splice!(P, i, filter(!isempty, [A, [x], N]))
                S, L = smaller_larger(A, N)
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

node_count(t::StrongModuleTree) = length(t)
node_count(x::Any) = 1

leaf_count(t::StrongModuleTree) = sum(leaf_count, t.nodes)
leaf_count(v::Vector) = sum(leaf_count, v)
leaf_count(x::Any) = 1

first_leaf(t::StrongModuleTree) = first_leaf(first(t.nodes))
first_leaf(v::Vector) = first_leaf(first(v))
first_leaf(x::Any) = x

last_leaf(t::StrongModuleTree) = last_leaf(last(t.nodes))
last_leaf(v::Vector) = last_leaf(last(v))
last_leaf(x::Any) = x

rand_leaf(t::StrongModuleTree) = rand_leaf(rand(t.nodes))
rand_leaf(v::Vector) = rand_leaf(last(v))
rand_leaf(x::Any) = x

function leaves(t::StrongModuleTree{T}) where T
    L = T[]
    for x in t.nodes
        append!(L, leaves(x))
    end
    return L
end
leaves(x::Any) = [x]

edge_string(t::StrongModuleTree, post::String="") =
    edge = t.kind == :prime    ? "" :
           t.kind == :complete ? "$(t.edge[1])$post" :
                                 "$(join(t.edge,"/"))$post"

kind_string(t::StrongModuleTree) = "$(edge_string(t,"-"))$(t.kind)"
Base.summary(t::StrongModuleTree) = "$(length(t))-node $(kind_string(t)) $(typeof(t))"

function Base.show(io::IO, t::StrongModuleTree)
    if get(io, :compact, false)
        print(io,
            edge_string(t,"-"), t.kind, " ",
            node_count(t), "-node (",
            leaf_count(t), "-leaf) module: ",
            first_leaf(t)
        )
    else
        parens = t.kind == :prime ? "{}" : t.kind == :linear ? "[]" : "()"
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

function Base.sort!(t::StrongModuleTree; lt=isless, by=first_leaf, rev::Bool=false)
    for x in t.nodes
        x isa StrongModuleTree || continue
        sort!(x, lt=lt, by=by, rev=rev)
    end
    sort!(t.nodes, lt=lt, by=by, rev=rev)
    return t
end

function StrongModuleTree(
    G::AbstractMatrix,
    v::AbstractVector{T},
    op::Vector{Int},
    cl::Vector{Int},
) where T
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
            t.kind == x.kind != :prime && t.edge == x.edge || continue
            splice!(t.nodes, i, x.nodes)
            i += length(x)
        end
    end

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
    return t
end

function StrongModuleTree(G::AbstractMatrix, p::Vector{Int})
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
for _ = 1:100
    # global n, G, p, T, p′, T′
    n = rand(3:10)
    G = rand(0:1, n, n)
    G .= G .⊻ G'
    @assert G == G'
    p = graph_factorizing_permutation(G)
    @assert is_modular_permutation(G, p)
    T = sort!(StrongModuleTree(G, p))
    @assert sort!(strong_modules(T), lt=lexless) ==
            sort!(strong_modules(G), lt=lexless)
    for _ = 1:10
        p′ = graph_factorizing_permutation(G, shuffle(1:n))
        @assert is_modular_permutation(G, p′)
        T′ = sort!(StrongModuleTree(G, p′))
        @assert T == T′
    end
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
            emit(x:y)
        end
    end
end

function all_common_intervals(p::Vector{Int})
    intervals = UnitRange{Int}[]
    all_common_intervals(p) do r
        push!(intervals, r)
    end
    return intervals
end

function strong_common_intervals(p::Vector{Int})
    intervals = all_common_intervals(p)
    return filter(A -> all(B -> !overlap(A, B), intervals), intervals)
end

permutation_graph(p1::Vector{Int}, p2::Vector{Int}) =
    Int[xor(p1[i] < p1[j], p2[i] < p2[j]) for i=1:length(p1), j=1:length(p2)]

## McConnell & Montgolfier 2004: "Linear-time modular decomposition of directed graphs"

## perfect tournament factorization ##

is_tournament(G::AbstractMatrix) = G + G' + I == ones(G)

function tournament_factorizing_permutation(
        T::AbstractMatrix,
        V::Vector{Int} = collect(1:checksquare(T)),
    )
    n, P = length(V), [V]
    for x = V
        i = findfirst(C->x in C, P)
        C = P[i]
        B = filter(y->x != y && T[x,y] < T[y,x], C)
        A = filter(y->x != y && T[y,x] < T[x,y], C)
        splice!(P, i, filter(!isempty, [B, [x], A]))
    end
    return map(first, P)
end

false &&
for _ = 1:1000
    # global n, T, p
    n = rand(3:10)
    T = Int[i != j && rand(Bool) for i=1:n, j=1:n]
    T .= T .⊻ T' .⊻ tril(ones(Int,n,n),-1)
    @assert is_tournament(T)
    p = tournament_factorizing_permutation(T)
    modules = all_modules(T)
    @assert is_modular_permutation(G, p, modules=modules)
end

## digraph factorization ##

function nodes!(v::Vector{StrongModuleTree{T}}, t::StrongModuleTree{T}) where T
    for x in t.nodes
        x isa StrongModuleTree || continue
        nodes!(v, x)
    end
    push!(v, t)
    return v
end
nodes(t::StrongModuleTree) = nodes!(typeof(t)[], t)

strong_modules(t::StrongModuleTree) = map(sort!∘leaves, nodes(t))

function parent_node(t::StrongModuleTree, S::Vector)
    for x in t.nodes
        x isa StrongModuleTree || continue
        L = leaves(x)
        isempty(S ∩ L) && continue
        S ⊆ L && length(S) < length(L) && return parent_node(x, S)
        break
    end
    return t
end

# TODO: implement linear time algorithm from this thesis:
## Dalhaus 1998: "Parallel algorithms for hierarchical clustering and
##   applications to split decomposition and parity graph recognition"
## There's also a 2000 paper by the same name, but the PDF is jumbled

function overlap_components(
    s::StrongModuleTree,
    t::StrongModuleTree,
    M = strong_modules(s),
    N = strong_modules(t),
)
    O = M ∪ N
    n = length(O)
    R = speye(Int, n)
    for i = 1:n-1, j = i+1:n
        overlap(O[i], O[j]) || continue
        R[i,j] = R[j,i] = 1
    end
    while true
        R′ = min.(1, R^2)
        R′ == R && break
        R .= R′
    end
    for (i, j) in zip(findn(R)...)
        i < j || continue
        O[i] = O[j] = sort!(O[i] ∪ O[j])
    end
    return unique(O)
end

function intersect_permutation(
    V::AbstractVector{E},
    s::StrongModuleTree{E},
    t::StrongModuleTree{E},
) where E
    Ms = strong_modules(s)
    Mt = strong_modules(t)
    U = filter!(overlap_components(s, t, Ms, Mt)) do X
        (X in Ms || parent_node(s, X).kind != :prime) &&
        (X in Mt || parent_node(t, X).kind != :prime)
    end
    for x in V; push!(U, [x]); end
    R = Dict()
    for X in U
        S = parent_node(t, X)
        T = parent_node(s, X)
        union!(get!(()->Set{Int}(), R, (S, T)), X)
    end
    N = U ∪ map(sort!∘collect, values(R))
    N = filter!(X->length(X) > 1, N)
    T = Any[[] for x in V]
    for node in sort!(N, by=length, rev=true)
        an = Vector{Any}(node)
        for x in node
            push!(T[x], an)
        end
    end
    for x in V, i = 1:length(T[x])-1
        parent = T[x][i]
        child = T[x][i+1]
        child in parent && continue
        filter!(y->!(y in child), parent)
        push!(parent, child)
    end
    T = T[1][1]
    p = E[]
    record_vals(v::Vector) = foreach(record_vals, v)
    record_vals(x::E) = push!(p, x)
    record_vals(T)
    return p
end

intersect_permutation(n::Int, s::StrongModuleTree{Int}, t::StrongModuleTree{Int}) =
    intersect_permutation(1:n, s, t)

function digraph_factorizing_permutation(G::AbstractMatrix)
    Gs = G .| G'
    Gd = G .& G'
    H = Gs + Gd
    n = checksquare(G)
    ps = graph_factorizing_permutation(Gs)
    pd = graph_factorizing_permutation(Gd)
    s = StrongModuleTree(Gs, ps)
    t = StrongModuleTree(Gd, pd)
    p = intersect_permutation(n, s, t)
    h = StrongModuleTree(H, p)
    function sort_leaves!(h)
        for x in h.nodes
            x isa StrongModuleTree && sort_leaves!(x)
        end
        h.kind == :complete || return
        if h.edge == (1,1) # tournament node
            X = map(first_leaf, h.nodes)
            q = tournament_factorizing_permutation(G, X)
            o = Dict(x => i for (i, x) in enumerate(q))
            sort!(h.nodes, by=x->o[first_leaf(x)])
        else # 0/2-complete node
            d = Dict()
            c = (1:n)\leaves(h)
            for x in h.nodes
                i = first_leaf(x)
                push!(get!(d, G[c,i], []), x)
            end
            h.nodes .= [x for v in values(d) for x in v]
        end
    end
    sort_leaves!(h)
    return leaves(h)
end

false &&
for i = 1:1000
    # global n, G, p
    n = rand(3:10)
    G = Int[i != j && rand(Bool) for i=1:n, j=1:n]
    p = digraph_factorizing_permutation(G)
    @assert is_modular_permutation(G, p)
end

nothing
