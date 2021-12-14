using LightGraphs, SimpleWeightedGraphs, DataStructures

function limited_dijkstra_shortest_paths(g::AbstractGraph,
    srcs::Vector{U},
    distmx::AbstractMatrix{T}=weights(g);
    allpaths=false,
    trackvertices=false,
    limit::T=zero(T)
    ) where T <: Real where U <: Integer

    nvg = nv(g)
    dists = fill(typemax(T), nvg)
    parents = zeros(U, nvg)
    visited = zeros(Bool, nvg)

    pathcounts = zeros(UInt64, nvg)
    preds = fill(Vector{U}(), nvg)
    H = PriorityQueue{U,T}()
    # fill creates only one array.
    
    dolimit = (limit > zero(T))

    for src in srcs
        dists[src] = zero(T)
        visited[src] = true
        pathcounts[src] = 1
        H[src] = zero(T)
    end

    closest_vertices = Vector{U}()  # Maintains vertices in order of distances from source
    sizehint!(closest_vertices, nvg)

    while !isempty(H)
        u = dequeue!(H)

        if trackvertices
            push!(closest_vertices, u)
        end

        d = dists[u] # Cannot be typemax if `u` is in the queue
        
        if dolimit
            if d > limit
                break
            end
        end
        
        for v in outneighbors(g, u)
            alt = d + distmx[u, v]

            if !visited[v]
                visited[v] = true
                dists[v] = alt
                parents[v] = u

                pathcounts[v] += pathcounts[u]
                if allpaths
                    preds[v] = [u;]
                end
                H[v] = alt
            elseif alt < dists[v]
                dists[v] = alt
                parents[v] = u
                #615
                pathcounts[v] = pathcounts[u]
                if allpaths
                    resize!(preds[v], 1)
                    preds[v][1] = u
                end
                H[v] = alt
            elseif alt == dists[v]
                pathcounts[v] += pathcounts[u]
                if allpaths
                    push!(preds[v], u)
                end
            end
        end
    end

    if trackvertices
        for s in vertices(g)
            if !visited[s]
                push!(closest_vertices, s)
            end
        end
    end

    for src in srcs
        pathcounts[src] = 1
        parents[src] = 0
        empty!(preds[src])
    end

    dijkstrastate = LightGraphs.DijkstraState{T,U}(parents, dists, preds, pathcounts, closest_vertices)
    trim_dijkstra_state!(dijkstrastate, limit, skipinf=true)
    return dijkstrastate
end

# single vertice
limited_dijkstra_shortest_paths(g::AbstractGraph, src::Integer, distmx::AbstractMatrix=weights(g); allpaths=false, trackvertices=false, limit=0.) =
limited_dijkstra_shortest_paths(g, [src;], distmx; allpaths=allpaths, trackvertices=trackvertices, limit=limit)

function trim_dijkstra_state(dijkstra_state::LightGraphs.DijkstraState{T,U}, limit::T=zero(T); skipinf=true) where T <: Real where U <: Integer
    s = LightGraphs.DijkstraState(
        deepcopy(dijkstra_state.parents),
        deepcopy(dijkstra_state.dists),
        deepcopy(dijkstra_state.predecessors), # no need for deep copy because we are not modifying content of predecessors itself
        deepcopy(dijkstra_state.pathcounts),
        deepcopy(dijkstra_state.closest_vertices)
    )
    return trim_dijkstra_state!(s, limit; skipinf=true)
end

function trim_dijkstra_state!(dijkstra_state::LightGraphs.DijkstraState{T,U}, limit::T=zero(T); skipinf=true) where T <: Real where U <: Integer
    s = dijkstra_state
    dolimit = (limit > zero(T))
    
    if (dolimit)
        for i in 1:length(s.dists)
            d = s.dists[i]
            if d > limit && (d < typemax(T) || !skipinf)
                s.parents[i] = zero(U)
                s.dists[i] = typemax(T)
                s.predecessors[i] = Vector{U}()
                s.pathcounts[i] = zero(U)
            end
        end
    end
    return s
end