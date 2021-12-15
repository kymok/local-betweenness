using SparseArrays, LinearAlgebra

include("./limited_dijkstra.jl")

function sortperm_unique_zerofirst(values)
    p = sortperm(values, rev=true)
    sorted = values[p]

    unique_vals = []
    unique_perm = []
    ctr = 0
    for i in 1:length(sorted)
        if i == 1 || sorted[i] != sorted[i-1] # new value
            push!(unique_perm, [])
            ctr += 1
            push!(unique_perm[ctr], p[i])
            push!(unique_vals, sorted[i])
        else sorted[i] == sorted[i-1] # existing value
            push!(unique_perm[ctr], p[i])       
        end
    end
    
    unique_vals_zerofirst = []
    unique_perm_zerofirst = []
    
    for i in 1:length(unique_vals)
        v = unique_vals[i]
        if v > 0
            push!(unique_vals_zerofirst, unique_vals[i])
            push!(unique_perm_zerofirst, unique_perm[i])
        elseif v==0
            pushfirst!(unique_vals_zerofirst, unique_vals[i])
            pushfirst!(unique_perm_zerofirst, unique_perm[i])
        end
    end
    
    return unique_vals_zerofirst, unique_perm_zerofirst
end

function _limited_edge_betweenness_reducer(
        a::Tuple{Array{Array{Float64,1},1}, Array{SparseMatrixCSC{Float64,Int64},1}},
        b::Tuple{Array{Array{Float64,1},1}, Array{SparseMatrixCSC{Float64,Int64},1}}
    )::Tuple{Array{Array{Float64,1},1}, Array{SparseMatrixCSC{Float64,Int64},1}}
    return (a[1]+b[1], a[2]+b[2])
end

function limited_edge_betweenness_centrality(g::AbstractGraph,
    vs::AbstractVector=vertices(g),
    distmx::AbstractMatrix{T}=weights(g);
    normalize=true,
    endpoints=false,
    limits=[zero(T)]) where T <: Real

    n_v = nv(g)
    k = length(vs)
    isdir = is_directed(g)

    betweenness = zeros(n_v)
    edge_betweenness = spzeros(n_v, n_v)
    
    #sort
    (limits_sorted, limits_perm) = sortperm_unique_zerofirst(limits)
    
    (betweenness, edge_betweenness) = @distributed _limited_edge_betweenness_reducer for s in vs
    #for s in vs
        temp_betweenness = [zeros(n_v) for i in 1:length(limits_sorted)]
        temp_edge_betweenness = [spzeros(n_v, n_v) for i in 1:length(limits_sorted)]
        
        if degree(g, s) > 0  # this might be 1?
            state = limited_dijkstra_shortest_paths(g, s, distmx; allpaths=true, trackvertices=true, limit=limits_sorted[1])
            
            for i in 1:length(limits_sorted)
                trim_dijkstra_state!(state, limits_sorted[i])
                if endpoints
                    _edge_accumulate_endpoints!(temp_betweenness[i], temp_edge_betweenness[i], state, g, s)
                else
                    _edge_accumulate_basic!(temp_betweenness[i], temp_edge_betweenness[i], state, g, s)
                end
            end
            
        end
        (temp_betweenness, temp_edge_betweenness)
    end
    
    if !isdir
        edge_betweenness = [SparseMatrixCSC(UpperTriangular(i + i')) for i in edge_betweenness]
    end
    
    # rescale
    for i in 1:length(betweenness)
        LightGraphs._rescale!(betweenness[i],
        n_v,
        normalize,
        isdir,
        k)
    end
    
    for i in 1:length(edge_betweenness)
        _edge_rescale!(edge_betweenness[i],
        n_v,
        normalize,
        isdir,
        k)
    end
    
    #unsort    
    betweenness_unsorted::Vector{Any} = [nothing for limit in limits]
    edge_betweenness_unsorted::Vector{Any} = [nothing for limit in limits]
    for (p, i) in zip(limits_perm, 1:length(limits_perm))
        for perm in p
            betweenness_unsorted[perm] = betweenness[i]
            edge_betweenness_unsorted[perm] = edge_betweenness[i]
        end
    end
    
    return (betweenness_unsorted, edge_betweenness_unsorted)
end


function limited_betweenness_centrality(g::AbstractGraph,
    vs::AbstractVector=vertices(g),
    distmx::AbstractMatrix{T}=weights(g);
    normalize=true,
    endpoints=false,
    limits=[zero(T)]) where T <: Real

    n_v = nv(g)
    k = length(vs)
    isdir = is_directed(g)

    betweenness = zeros(n_v)
    edge_betweenness = spzeros(n_v, n_v)
    
    #sort
    (limits_sorted, limits_perm) = sortperm_unique_zerofirst(limits)
    
    betweenness = @distributed (+) for s in vs
    #for s in vs
        temp_betweenness = [zeros(n_v) for i in 1:length(limits_sorted)]
        
        if degree(g, s) > 0  # this might be 1?
            state = limited_dijkstra_shortest_paths(g, s, distmx; allpaths=true, trackvertices=true, limit=limits_sorted[1])
            
            for i in 1:length(limits_sorted)
                trim_dijkstra_state!(state, limits_sorted[i])
                if endpoints
                    LightGraphs._accumulate_endpoints!(temp_betweenness[i], state, g, s)
                else
                    LightGraphs._accumulate_basic!(temp_betweenness[i], state, g, s)
                end
            end
        end
        temp_betweenness
    end
        
    # rescale
    for i in 1:length(betweenness)
        LightGraphs._rescale!(betweenness[i],
        n_v,
        normalize,
        isdir,
        k)
    end
    
    #unsort    
    betweenness_unsorted::Vector{Any} = [nothing for limit in limits]
    for (p, i) in zip(limits_perm, 1:length(limits_perm))
        for perm in p
            betweenness_unsorted[perm] = betweenness[i]
        end
    end
    
    return betweenness_unsorted
end


# Accumulation

function _edge_accumulate_basic!(betweenness::Vector{Float64},
    edge_betweenness::AbstractMatrix{Float64},
    state::LightGraphs.DijkstraState,
    g::AbstractGraph,
    si::Integer)

    n_v = length(state.parents) # this is the ttl number of vertices
    δ = zeros(n_v)
    σ = state.pathcounts
    P = state.predecessors

    # make sure the source index has no parents.
    P[si] = []
    # we need to order the source vertices by decreasing distance for this to work.
    S = reverse(state.closest_vertices) #Replaced sortperm with this
    for w in S
        coeff = (1.0 + δ[w]) / σ[w]
        for v in P[w]
            if v > 0
                δ[v] += (σ[v] * coeff)
                edge_betweenness[v, w] += (σ[v] * coeff)
            end
        end
        if w != si
            betweenness[w] += δ[w]
        end
    end
    return nothing
end

function _edge_accumulate_endpoints!(betweenness::Vector{Float64},
    edge_betweenness::AbstractMatrix{Float64},
    state::LightGraphs.DijkstraState,
    g::AbstractGraph,
    si::Integer)

    n_v = length(state.parents) # this is the ttl number of vertices
    δ = zeros(n_v)
    σ = state.pathcounts
    P = state.predecessors

    # make sure the source index has no parents.
    P[si] = []
    # we need to order the source vertices by decreasing distance for this to work.
    S = reverse(state.closest_vertices) #Replaced sortperm with this
    
    betweenness[s] += length(S) - 1    # 289
    
    for w in S
        coeff = (1.0 + δ[w]) / σ[w]
        for v in P[w]
            δ[v] += σ[v] * coeff
            edge_betweenness[v, w] += (σ[v] * coeff)
        end
        if w != si
            betweenness[w] += (δ[w] + 1)
        end
    end
    return nothing
    
end


# rescale

function _edge_rescale!(edge_betweenness::AbstractMatrix{Float64}, n::Integer, normalize::Bool, directed::Bool, k::Integer)
    if normalize
        if n <= 2
            do_scale = false
        else
            do_scale = true
            scale = 1.0 / (n * (n - 1))
        end
    else
        if !directed
            do_scale = true
            scale = 1.0 / 2.0
        else
            do_scale = false
        end
    end
    if do_scale
        if k > 0
            scale = scale * n / k
        end
        edge_betweenness .*= scale

    end
    return nothing
end