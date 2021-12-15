# local-betweenness

## About
This is the code for https://doi.org/10.1007/978-3-030-76059-5_26

## Usage

### Exporting a graph from Python

see `graph_tool_to_graphml.py`. Property name must be `e_weight` (edge weight) and `eid` (edge id).

### Computation
```bash
julia betweenness/src/betweenness.jl --input ./graphml/sample.graphml --output ./
```

The julia code takes a GraphML data exported from [`graph-tools`](https://graph-tool.skewed.de) as input. The code relies on the XML structure that the library exports (FIXME).

## License

### Code
- MIT otherwise noted.

### Road Network Data
- Data © OpenStreetMap contributors.
- The graph data was retrieved using [OSMnx](https://github.com/gboeing/osmnx) from the OpenStreetMap.