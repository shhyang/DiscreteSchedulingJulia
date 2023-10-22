module DiscreteScheduling
using Graphs,LinearAlgebra,SparseArrays
using CDDLib,Polyhedra
include("./MaxIndependentSets.jl")

using .MaxIndependentSets

export WeightedGraph, MaxSchedulingGraph, SchedulingGraph
export maxindepset, linenet, extendgraph, schedgraph, maxschedgraph
export RateRegion, genVF, OptimalRate

"""
A network model with a directed graph G and a delay matrix D. 
L is the number of links
"""
struct WeightedGraph
    L::Int
    G::SimpleDiGraph
    D::Matrix{Int64}
end

"""
    linenet(L,K)

generate a line network of length L with the K-hop interference model
"""
function linenet(L,K)
    G = SimpleDiGraph(L)
    D = Matrix{Int64}(undef,L,L)
    for i = 1 : L 
        for j = 1 : L
            D[i,j] = 1 - abs(j-i-1)
            if (j != i) &&  ( abs(j-i-1) <= K )
                add_edge!(G,i,j)
            end    
        end
    end
    return WeightedGraph(L,G,D)
end

"""
    extendgraph(wg,T)

generate a graph by extending the weighted graph wg by T times
"""
function extendgraph(wg::WeightedGraph,T::Int)
    L = wg.L
    G = SimpleGraph(L*T) # (l,t) is L*(t-1)+l 
    a = (l,t) -> L*(t-1)+l
    for l = 1 : L
        for t = 1 : T
            for l2 in neighbors(wg.G,l)      
                t2 = t + wg.D[l,l2]
                if t2>=1 && t2<=T
                    add_edge!(G, a(l,t), a(l2,t2))
                end
            end
        end
    end
    return G
end

"""
    MaxSchedulingGraph

The max scheduling graph (ML*, MR*, E*)
"""
struct MaxSchedulingGraph
    T::Int
    ML::Array{Array}
    MR::Array{Array}
    A::SparseMatrixCSC{Bool,Int64}
end

"""
    maxschedgraph(wg,T)

Generate the max step-T scheduling graph of wg 
"""
function maxschedgraph(wg::WeightedGraph, T::Int)
    L = wg.L
    indepset = maxindepset(extendgraph(wg, T+T), 3)
    sML = Array[]
    sMR = Array[]
    ML = Matrix{Bool}[]
    MR = Matrix{Bool}[]
    A = spzeros(Bool,2^(L*T),2^(L*T))
    function add!(M,v)
        m = zeros(Bool,L,T)
        for b in v
            l = (b-1)%L+1
            t = ceil(Int,b/L)
            m[l,t] = 1
        end
        push!(M,m)
    end
    for E in indepset
        vl = sort(intersect(E,Set(1:L*T)))
        i = findfirst(x -> x == vl, sML)
        if isnothing(i)
            push!(sML, vl)
            i = length(sML)
            add!(ML,vl)
        end
        vr = sort(intersect(E,Set(L*T+1:L*(T+T)))) .- L*T
        j = findfirst(x -> x == vr, sMR)
        if isnothing(j)
            push!(sMR, vr)
            j = length(sMR)
            add!(MR,vr)
        end
        A[i,j] = 1
    end
    return MaxSchedulingGraph(T,ML,MR,A[1:length(ML),1:length(MR)])
end

"""
    SchedulingGraph

The scheduling graph (M_T, E_T,Q)
"""
struct SchedulingGraph
    T::Int
    Q::Int
    V::Array{Array}
    F::DiGraph
end

"""
    schedugraph(wg,T,Q)

Generate the step-Q, length-T scheduling graph of wg. 
See V-A in [1], Algorithm 1
"""
function schedgraph(wg::WeightedGraph, T::Int, Q::Int)
    L = wg.L
    indepset = maxindepset(extendgraph(wg, T+Q),3)
    V = Array[]
    F = DiGraph(2^(L*T))
    function addE!(E)
        vl = sort(intersect(E,Set(1:L*T)))
        vr = sort(intersect(E,Set(L*Q+1:L*(T+Q)))) .- L*Q
        i = findfirst(x -> x == vl, V)
        j = findfirst(x -> x == vr, V)
        if isnothing(i) || isnothing(j)
            if isnothing(i)
                push!(V,vl)
                i = length(V)
                j = findfirst(x -> x == vr, V)
            end      
            if isnothing(j)
                push!(V,vr)
                j = length(V)
            end
        elseif adjacency_matrix(F)[i,j] == 1
            return
        end
        add_edge!(F,i,j)
        for e in E  
            addE!(setdiff(E,Set([e])))
        end
    end
    for E in indepset
        addE!(E)
    end

    for i in eachindex(V)
        m = zeros(Bool,L,T)
        for b in V[i]
            l = (b-1)%L+1
            t = ceil(Int,b/L)
            m[l,t] = 1
        end
        V[i] = m
    end

    return SchedulingGraph(T,Q, V, F[1:length(V)])
end

"""
    genVF(sg)

Generate the graph (V, F) from the max scheduling graph sg.
See V-D in [1].
"""
function genVF(sg::MaxSchedulingGraph)
    V = Array[]
    F = DiGraph(length(sg.ML)*length(sg.MR))
    for B1 in sg.MR
        for A3 in sg.ML
            for a in findall(sg.A)
                v1 = B1 .& sg.ML[a[1]]
                v2 = sg.MR[a[2]] .& A3
                i = findfirst(x -> x == v1, V)
                if isnothing(i)
                    push!(V, v1)
                    i = length(V)
                end      
                j = findfirst(x -> x == v2, V)
                if isnothing(j)
                    push!(V, v2)
                    j = length(V)
                end
                add_edge!(F,i,j)
            end
        end
    end
    return SchedulingGraph(sg.T,sg.T, V,F[1:length(V)])
end

function pushtoR(S,r)
    for b in S
        if iszero(.!(r .<= b))
            return
        end
        if iszero(.!(b .<= r))
            delete!(S,b)
        end
    end    
    push!(S,r)
end

"""
    RateRegion(sg,k)

Calculate the rate region formed by cycles of length up to K in the scheduling graph sg.
The same algorithm can apply on (V, F). 
"""
function RateRegion(sg::SchedulingGraph,K)
    cycles = simplecycles_limited_length(sg.F,K)
    RR = Set{Vector{Float64}}()
    for cyc in cycles
        v = sum(sum(sg.V[cyc]),dims=2)[:]./length(cyc)
        pushtoR(RR, v)
    end
    r = collect(RR)./sg.T
    pp = polyhedron(vrep(r),CDDLib.Library())
    removevredundancy!(pp)
    return pp
end

"""
    W2AB(sg)

Call calculation of W2(A,B) for all A in ML* and B in MR* from the max scheduling graph sg
See V-C in [1], Algorithm 2
"""
function W2AB(sg::MaxSchedulingGraph)
    l1 = length(sg.ML)
    l2 = length(sg.MR)
    R2 = Array{Set}(undef,l1,l2)
    for i = 1 : l1
        for j = 1 : l2
            R2[i,j] = Set{Vector{Int64}}()
            for i2 in findall(!iszero, sg.A[:,j])
                for j2 in findall(!iszero, sg.A[i,:])
                    B1A2 = sum(sg.MR[j2] .& sg.ML[i2], dims=2)[:]
                    pushtoR(R2[i,j], B1A2)
                end
            end
        end
    end
    return R2
end


"""
    WAB(R, sg)

Call calculation of Wk(A,B) for all A in ML* and B in MR* with input Wk-1 and the max scheduling graph sg 
See V-C in [1], Algorithm 2
"""
function WAB(R,sg::MaxSchedulingGraph)
    l1 = length(sg.ML)
    l2 = length(sg.MR)
    R2 = Array{Set}(undef,l1,l2)
    for i = 1 : l1
        for j = 1 : l2
            R2[i,j] = Set{Vector{Int64}}()
            for i2 in findall(!iszero, sg.A[:,j])
                for j2 = 1 : l2
                    B1A2 = sum(sg.MR[j2] .& sg.ML[i2], dims=2)[:]
                    if isempty(R[i,j2])
                        pushtoR(R2[i,j], b)
                    else
                        for b in R[i,j2]
                            pushtoR(R2[i,j], B1A2+b)
                        end
                    end
                end
            end
        end
    end
    return R2
end

"""
    RateRegion(sg,K)

Calculate the region region using the max scheduling graph sg up to cycles of length K.
See Algorithm 3 in [1]
"""
function RateRegion(sg::MaxSchedulingGraph,K)
    #RR = Array{Array}(undef,K)
    ll = length(sg.ML)
    lr = length(sg.MR)
    RR = Set{Vector{Float64}}()
    for i = 1 : length(sg.ML)
        for j in findall(!iszero, sg.A[i,:])
            pushtoR(RR,sum(sg.ML[i] .& sg.MR[j], dims=2)[:]./sg.T)
        end
    end
    #RR[1] = collect(R1)./sg.T
    function calRR(R,k)
        for r = 1 : lr
            for l = 1 : ll
                AB = sum(sg.MR[r] .& sg.ML[l], dims=2)[:]
                for b in R[l,r]
                    pushtoR(RR, (AB + b)./(k*sg.T))
                end
            end
        end
    end
    if K > 1
        R2 = W2AB(sg)
        calRR(R2,2)
        for k=3:K
            Rplus = WAB(R2,sg)
            calRR(Rplus,k)
            R2 = Rplus
        end
    end
    #return RR     
    p = polyhedron(vrep(collect(RR)),CDDLib.Library())
    removevredundancy!(p)
    return p
end

"""
    U2AB(sg,f)

Calculate U2 for max scheduling graph sg and function f.
See Algorithm 4 in [1]
"""
function U2AB(sg::MaxSchedulingGraph,f::Function)
    l1 = length(sg.ML)
    l2 = length(sg.MR)
    U2 = zeros(Float64,l1,l2)
    for i = 1 : l1
        for j = 1 : l2
            #C2 = zeros(Bool,size(sg.ML[1]))
            for i2 in findall(!iszero, sg.A[:,j])
                for j2 in findall(!iszero, sg.A[i,:])
                    r = f(sum(sg.MR[j2] .& sg.ML[i2],dims=2))
                    if U2[i,j] < r
                        U2[i,j] = r
                        #C2 = sg.MR[j2] .& sg.ML[i2]
                    end
                end
            end
        end
    end
    return U2
end

"""
    UAB(sg,f)

Calculate Uk from Uk-1 for max scheduling graph sg and function f.
See Algorithm 4 in [1]
"""
function UAB(U,sg::MaxSchedulingGraph,f::Function)
    l1 = length(sg.ML)
    l2 = length(sg.MR)
    U2 = zeros(Float64,l1,l2)
    for i = 1 : l1
        for j = 1 : l2
            #C2 = zeros(Bool,size(sg.ML[1]))
            for i2 in findall(!iszero, sg.A[:,j])
                for j2 = 1 : l2
                    r = U[i,j2] + f(sum(sg.MR[j2] .& sg.ML[i2],dims=2))
                    if U2[i,j] < r
                        U2[i,j] = r
                        #C2 = sg.MR[j2] .& sg.ML[i2]
                    end
                end
            end
            #push!(C[i,j],C2)
        end
    end
    return U2
end

"""
    OptimalRate(f,sg,K)

Calcuate rate the optimal value of function f over the rate region up to cycle length K.
The function also returns a cycle that achieves the optimal value. 
See Algorithm 5 and 6 in [1]
"""
function OptimalRate(f::Function,sg::MaxSchedulingGraph,K)
    ll = length(sg.ML)
    lr = length(sg.MR)
    U = Array{Matrix}(undef,K)
    UU = Array{Float64}(undef,K)
    UU[1] = 0.0
    for i = 1 : ll
        for j in findall(!iszero, sg.A[i,:])
            u = f(sum(sg.ML[i] .& sg.MR[j],dims=2))
            if UU[1] < u
                UU[1] = u
            end
        end
    end
    UU[1] = UU[1]/sg.T

    function maxU(uk)
        ru = 0.0
        for i = 1 : ll
            for j = 1 : lr
                r = uk[i,j] + f(sum(sg.ML[i] .& sg.MR[j],dims=2))
                if ru < r
                    ru = r
                end
            end
        end
        return ru
    end

    U[2] = U2AB(sg,f)
    UU[2] = maxU(U[2])/2/sg.T

    for k=3:K
        U[k] = UAB(U[k-1],sg,f)
        UU[k] = maxU(U[k])/k/sg.T
    end

    # find the optimal cycle
    a = Array{Int}(undef,K)
    b = Array{Int}(undef,K)
    C = Array{Matrix}(undef,K)
    for i = 1 : ll
        for j = 1 : lr
            uu = U[K][i,j] + f(sum(sg.ML[i] .& sg.MR[j],dims=2))
            if (uu == UU[K]*K*sg.T)
                a[1] = i
                b[K] = j
                C[K] = sg.ML[i] .& sg.MR[j]
                break
            end
        end
    end

    for k = K:-1:3
        for i in findall(!iszero, sg.A[:,b[k]])
            for j = 1 : lr
                uu = U[k-1][a[1],j] + f(sum(sg.ML[i] .& sg.MR[j],dims=2))
                if (uu == U[k][a[1],b[k]])
                    a[k] = i
                    b[k-1] = j
                    C[k-1] = sg.ML[i] .& sg.MR[j]
                    break
                end
            end
        end
    end

    for i in findall(!iszero, sg.A[:,b[2]])
        for j in findall(!iszero, sg.A[a[1],:])
            if (U[2][a[1],b[2]] == f(sum(sg.ML[i] .& sg.MR[j],dims=2)))
                a[2] = i
                b[1] = j
                C[1] = sg.ML[i] .& sg.MR[j]
                break
            end
        end
    end
    return UU, C
end 

end