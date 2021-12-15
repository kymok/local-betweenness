using ArgParse
using LightGraphs, SimpleWeightedGraphs, GraphIO, EzXML
using Distributed
using DataFrames, CSV

LIMITS = collect(300:100:5000)

include(joinpath(@__DIR__, "jl/parse_graphml.jl"))

function parse_commandline()
    s = ArgParseSettings()
    @add_arg_table! s begin
        "--input"
            help = "input file path"
            required = true
        "--output"
            help = "output file path or directory"
            required = true
    end
    return parse_args(s)
end

function _splitdir(path)
    if(isdir(path))
        return (path, "")
    end
    return splitdir(path)
end

parsed_args = parse_commandline()
input_dir, input_file = _splitdir(abspath(parsed_args["input"]))
output_dir, output_file = _splitdir(abspath(parsed_args["output"]))
output_file = output_file == "" ? input_file : output_file
output_file, _ = splitext(output_file)
println("Input:  " * input_dir * " " * input_file)
println("Output: " * output_dir * " " * output_file)

doc = EzXML.readxml(open(joinpath(input_dir, input_file)))
g = parse_graphtool_GraphML(doc)

# Compute Betweenness
_procs = addprocs();
_g = g[:graph];
limits = convert(Array{Float64, 1}, LIMITS)
@everywhere include(joinpath(@__DIR__, "jl/distributed_limited_edge_betweenness.jl"))
@time begin
    betweenness, edge_betweenness = limited_edge_betweenness_centrality(_g, normalize=false, limits=limits)
end
rmprocs(_procs)


# Write Data
# Write Edge Betweenness 
eb_coords = reduce(unique∘vcat, [findall(!iszero, eb) for eb in edge_betweenness])
eb_data = reduce(hcat, [[edge_betweenness[i][c] for c in eb_coords] for i in 1:length(edge_betweenness)])

df_e = DataFrame(eb_data, :auto)
rename!(df_e, map(Symbol∘string, limits))
df_e.src_vid = [c[1]-1 for c in eb_coords]
df_e.dst_vid = [c[2]-1 for c in eb_coords]
output_file_e = output_file * "_edge.csv"
CSV.write(joinpath(output_dir, output_file_e), df_e)


# Write Vertex Betweenness
df_v = DataFrame(betweenness, :auto)
rename!(df_v, map(Symbol∘string, limits))
df_v.vid = [i-1 for i in 1:size(df_v)[1]]
output_file_v = output_file * "_vertex.csv"
CSV.write(joinpath(output_dir, output_file_v), df_v)
