# %%
from graph_tool.all import *

# %%
# Generate a small weighted graph
g = Graph(directed=False)
g.ep["e_weight"] = g.new_edge_property("double")
g.ep["eid"] = g.new_edge_property("int")
n = 5
g.add_vertex(n)
for i in range(n-1):
    e = g.add_edge(g.vertex(i), g.vertex(i+1))
    g.ep.e_weight[e] = i + 1
    g.ep.eid[e] = i
# %%
# export to a GraphML file
g.save("sample.graphml")

# %%
