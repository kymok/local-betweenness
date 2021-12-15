using LightGraphs, GraphIO
using SimpleWeightedGraphs
using EzXML

function parse_graphtool_GraphML(doc)
    
    g = SimpleWeightedGraph(SimpleGraph(0))
    
    weight_name = "e_weight"
    weight_key_id = ""
    eid_name = "eid"
    eid_key_id = ""
    nodeid2index = Dict{String, Int64}()
    edgeid2index = Dict{String, Int64}()
    eid2index = Dict{String, Int64}()
    
    # TODO can be improved using:
    # parse.nodeids="canonical" parse.edgeids="canonical" parse.order="nodesfirst"
    for i in eachelement(root(doc))
        
        if (i.name == "key")
            key_id = i["id"]
            key_name = i["attr.name"]
            if (key_name == weight_name)
                weight_key_id = key_id
            elseif (key_name == eid_name)
                eid_key_id = key_id
            end
        end
        
        if (i.name == "graph")
            #node
            for elem in eachelement(i)
                
                # First pass: add nodes
                if (elem.name == "node")
                    n_nodes = length(nodeid2index)
                    node_id = elem["id"]
                    nodeid2index[node_id] = n_nodes + 1
                    add_vertex!(g)
                end
            end
            
            #edge
            for elem in eachelement(i)
                
                # Second pass: add edges
                if (elem.name == "edge")
                    src_id = nodeid2index[elem["source"]]
                    trg_id = nodeid2index[elem["target"]]
                    edge_id = elem["id"]

                    weight = 0.
                    eid = 0
                    has_weight = false
                    for data in eachelement(elem)
                        if (data.name == "data")
                            edata_key = data["key"]
                            if (edata_key == weight_key_id)                                
                                has_weight = true
                                weight = parse(Float64, data.content)
                            elseif (edata_key == eid_key_id)
                                eid = data.content
                            end
                        end
                    end
                    
                    if (has_weight)
                        n_edges = length(edgeid2index)
                        edgeid2index[edge_id] = n_edges + 1
                        eid2index[eid] = n_edges + 1
                        add_edge!(g, src_id, trg_id, weight)
                    end
                end
            end
        end
    end

    return Dict(:graph=>g, :node_id=>nodeid2index, :edge_id=>edgeid2index, :eid=> eid2index)
end