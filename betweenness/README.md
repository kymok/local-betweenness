# Local Betweenness

Modified code of original `JuliaGraphs`.

## Usage

### Set limits

`src/betweenness.jl` Line 8:

default

```julia
LIMITS=collect(300:100:5000)
```

no limit

```julia
LIMITS = [0]
```

### Run

Following command processes a graphml file `input.graphml` and saves the result to `output_dir/input_{edge,vertex}.csv`:

```bash
julia src/betweenness.jl --input (input.graphml) output (output_dir)
```

### Batch Run

```
src/compute_all.sh input_dir output_dir
```

This command processes all `.graphml` files in `input_dir` and writes output to `output_dir`.


## Known Issues

### Fails when RAM is insufficient

Currently the program fails silently when RAM is insufficient. If some results are missing, try increasing RAM or executing with fewer number of parallel processes. (see Julia's `addprocs()` for further information)

### Inefficient computation when larage number of `nprocs`

