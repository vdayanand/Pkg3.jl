## Some graph functions ##

#=
X = sparse([1,1,2,2,3,4,4], [2,5,5,3,4,5,6], 1, 6, 6)
X += X'
G = dropzeros!(1 - X)
=#

N(G, v) = find(G[:, v])
const \ = setdiff

function BronKerboschTomita(emit, G, R, P, X)
    let available = unique(pkg_map[v] for V in (R, P) for v in V)
        # if any in R are unsatisfiable, return
        for v in R, r in req_map[v]
            r in available || return
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
    isempty(P) && isempty(X) && (emit(R); return)
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
        BronKerboschTomita(emit, G, [R; v], P \ Nv, X \ Nv)
        filter!(x -> x != v, P)
        push!(X, v)
    end
end

function maximal_indepedent_sets(io::IO, G::AbstractMatrix, inds::Vector{Int} = collect(1:size(G,2)))
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
maximal_indepedent_sets(path::String, G::AbstractMatrix, inds::Vector{Int} = collect(1:size(G,2))) =
    open(io->maximal_indepedent_sets(io, G, inds), path, "w")
maximal_indepedent_sets(G::AbstractMatrix, inds::Vector{Int} = collect(1:size(G,2))) =
    maximal_indepedent_sets(STDOUT, G, inds)

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

overlap(A::Vector, B::Vector) = !isempty(A \ B) && !isempty(A ∩ B) && !isempty(B \ A)

## McConnell & Montgolfier 2004: "Linear-time modular decomposition of directed graphs"

#=
n = 5
G = rand(n, n)
T = float(G .< G')
@assert T + T' + I == ones(T)
p = tournament_factorizing_permutation(T)
modules = filter(S->is_module(T, S), collect(subsets(1:n)))
@assert all(M->all(x->x == 1, diff(findin(M, p))), modules)
=#

tournament_factorizing_permutation(T::AbstractMatrix) =
    tfp!(T, collect(1:Base.LinAlg.checksquare(T)))

function tfp!(T::AbstractMatrix, p::Vector{Int}, lo::Int=1, hi::Int=length(p))
    if hi - lo > 1
        v = p[lo]
        i, j = lo + 1, hi
        while true
            while i < j && T[v,p[i]] == 0; i += 1; end;
            while i < j && T[p[j],v] == 0; j -= 1; end;
            i < j || break
            p[i], p[j] = p[j], p[i]
        end
        p[j], p[lo] = v, p[j]
        tfp!(T, p, lo, j-1)
        tfp!(T, p, j+1, hi)
    end
    return p
end

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

#=
n = 6
G = Int[i != j && rand() < 0.5 for i = 1:n, j = 1:n]
G .= G .⊻ G'
@assert G == G'
p = graph_factorizing_permutation(G)
modules = filter!(S->length(S) > 1 && is_module(G, S), collect(subsets(1:length(p))))
@assert all(M->all(x->x == 1, diff(findin(M, p))), modules)
strong = filter(A -> all(B -> !overlap(A, B), modules), modules)
=#

graph_factorizing_permutation(G::AbstractMatrix) =
    gfp!(G, collect(1:Base.LinAlg.checksquare(G)))

function gfp!(G::AbstractMatrix, p::Vector{Int}, lo::Int=1, hi::Int=length(p))
    if hi - lo > 1
        x = p[lo]
        i, j = lo+1, hi
        while true
            while i < j && G[p[i],x] == 0; i += 1; end;
            while i < j && G[p[j],x] != 0; j -= 1; end;
            i < j || break
            p[i], p[j] = p[j], p[i]
            i += 1; j -= 1
        end
        gfp!(G, p, lo+1, j)
        gfp!(G, p, j+1, hi)
    end
    return p
end

## Capelle, Habib & de Montgolfier 2002: "Graph decompositions and factorizing permutations"

function make_tree(v::AbstractVector{Int}, op::Vector{Int}, cl::Vector{Int})
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
    return s[end]
end

function classify_nodes(G::AbstractMatrix, t::Vector)
    node(t::Vector) = node(first(t))
    node(x::Any) = x
    order(a, b) = a < b ? (a, b) : (b, a)
    # check consecutive edges
    prime = false
    x, y = node(t[1]), node(t[2])
    edge = order(G[y,x], G[x,y])
    for i = 2:length(t)-1
        x, y = node(t[i]), node(t[i+1])
        if edge != order(G[y,x], G[x,y])
            prime = true
            break
        end
    end
    (prime ? () : edge) => map(x->x isa Vector ? classify_nodes(G, x) : x, t)
end

function ftree(G::AbstractMatrix, p::Vector{Int}=graph_factorizing_permutation(G))
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
                l = i
                i < j || continue
                if i <= lc[j-1] < uc[j-1] <= k
                    if c > 0
                        op[i] += 1
                        cl[k] += 1
                        l = k + 1
                    end
                else
                    if i < j-1
                        op[i] += 1
                        cl[j-1] += 1
                    end
                    l = k + 1
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
    # construct the tree
    t = make_tree(p, op, cl)
    # classify_nodes!(G, t)
end

function print_tree(io::IO, pair::Pair)
    show(io, pair[1])
    print_tree(io, pair[2])
end
function print_tree(io::IO, tr::Vector)
    print(io, "[")
    for (i, x) in enumerate(tr)
        print_tree(io, x)
        i < length(tr) && print(io, " ")
    end
    print(io, "]")
end
print_tree(io::IO, x::Any) = show(io, x)
print_tree(x::Any) = print_tree(STDOUT, x)

## Uno & Yagiura 2000: "Fast algorithms to enumerate all common intervals of two permutations"

#=
Brute force all common intervals:
[(x,y) for x=1:n-1 for y=x+1:n for l=1:n-1 for u=l+1:n if sort!(p1[x:y]) == sort!(p2[l:u])]
=#

function cif(p, i, j)
    i < j || return -1
    l, u = extrema(p[k] for k = i:j)
    (u-l) - (j-i)
end

if false
    n = length(p)
    F = [cif(p, i, j) for i = 1:n, j = 1:n]
end

function all_common_intervals(emit::Function, p::Vector{Int})
    for x = 1:length(p)-1
        l = u = p[x]
        for y = x+1:length(p)
            v = p[y]
            l, u = min(v, l), max(v, u)
            y - x < u - l && continue
            y - x > u - l && break
            emit(x, y)
        end
    end
end

function all_common_intervals(p::Vector{Int})
    intervals = NTuple{2,Int}[]
    all_common_intervals(p) do x, y
        push!(intervals, (x, y))
    end
    return intervals
end

all_common_intervals(emit::Function, p1::Vector{Int}, p2::Vector{Int}) =
    all_common_intervals(emit, invperm(p2)[p1])
all_common_intervals(p1::Vector{Int}, p2::Vector{Int}) =
    all_common_intervals(invperm(p2)[p1])
