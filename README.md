# MetricGraphs

[![Build Status](https://github.com/mathmd/MetricGraphs.jl/actions/workflows/CI.yml/badge.svg?branch=main)](https://github.com/mathmd/MetricGraphs.jl/actions/workflows/CI.yml?query=branch%3Amain)
This objective of this Julia package is to provide numerical solvers for PDEs on metric graphs.
As a start, only first order accurate implicit finite discretization scheme is implemented for a non-homogeneous reaction diffusion equation with Neumann-Kirhhoff condition (See `rd_dynamics!` in `/src/MetricGraphs.jl`.)A usage example can be found in the folder `/examples/`.

This package also provide data type `EdgeVector` and `VertexVector` for easy indexing with `Edge` which is updated as the underlying graph structure is edited (Warning: adding more vertex or adding/removing edges do not cause a problem, but do not remove existing vertices.) 

WIP:

- Extending Neumann-Kirhhoff vertex conditions to $\delta$-type vertex conditions, including Dirichlet condition.
- Implementing higher order accurate implicit finite discretization scheme.

To use this code, we simply need to define a `MetricGraph` from a DiGraph (using `Graphs.jl`) and edge lengths as `EdgeVector`. 
Then, discretize it using `subdivide_graph` which inserts $N$ vertices on each existing edges; this will output a discretized `MetricGraph` subdivided lengths, along with `vertex_map` and `edge_map` which keep track of the discretized graph.
The discretized graph is used as the domain for the numerical method.

Future plan:

- Implementing pseudo spectral method
- Seamless integration with existing software for numerical continuation and bifurcation analysis.

## Installation

This package is not yet registered in the official Julia repository.
To add this to the local Julia packages, enter the `Pkg` mode with `]` in REPL and enter
```
dev git@github.com:mathmd/MetricGraphs.jl.git
add MetricGraphs
```
Then, `using MetricGraphs` normally to use this package (see examples in `/examples/`.)
