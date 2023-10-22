module MaxIndependentSets
using Graphs

export maxindepset

# Bron Kerbosch algorithm for maximal independent set

function BK1(pg,R,P,X,cliques)
    if isempty(P) && isempty(X)
        push!(cliques,collect(R)) 
        return 
    end
    for v in P
        N = neighbors(pg,v)
        BK1(pg,union(R,Set([v])),intersect(P,N),intersect(X,N),cliques)
        setdiff!(P,Set([v]))
        union!(X,Set([v]))
    end
end


function pivot(pg,P)
    p = collect(P)
    (v,idx) = findmin(sum(adjacency_matrix(pg)[p,p],dims=2))
    return p[idx]
end

function BK2(pg,R,P,X,cliques)
    if isempty(P) && isempty(X)
        push!(cliques,collect(R)) 
        return 
    end
    u = pivot(pg,union(P,X))
    for v in setdiff(P,neighbors(pg,u))
        N = neighbors(pg,v)
        BK2(pg,union(R,Set([v])),intersect(P,N),intersect(X,N),cliques)
        setdiff!(P,Set([v]))
        union!(X,Set([v]))
    end
end

function degeneracy(pg)
    l = nv(pg)
    dgc = maximum(core_number(pg))
    deg = Array(1:l)
    idx = Array(1:l)
    for i = 1 : l-1
        d = findfirst(x -> x <= dgc, degree(pg[idx]))
        deg[i] = idx[d]
        splice!(idx,d)
    end
    deg[l] = idx[1]
    return deg
end

function maxindepset(pg::SimpleGraph,t)
    g = complement(pg)
    indepset = Vector[] # 
    R = Set{Int}()
    X = Set{Int}()
    P = Set(1 : nv(pg)) 
    if t == 1
        BK1(g,R,P,X,indepset)
        return indepset 
    elseif t ==2
        BK2(g,R,P,X,indepset)
        return indepset
    end
    # t == 3
    for v in degeneracy(g)
        N = neighbors(g,v)
        BK2(g,Set([v]),intersect(P,N),intersect(X,N),indepset)
        setdiff!(P,Set([v]))
        union!(X,Set([v]))
    end
    return indepset 
end



end