module MetricGraphs

using Graphs

include("./Types.jl")
export Hamiltonian, DivergenceForm, MetricGraph

include("./Utils.jl")
export subdivide_graph

using LinearAlgebra

# export rd_dynamics@


end # module
