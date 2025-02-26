#%%
using MetricGraphs
using Graphs
using Plots

M = 5
G = DiGraph(4M + 5)
for i in 0:M-1
        add_edge!(G, 4i + 1, 4i + 2)
        add_edge!(G, 4i + 2, 4i + 3)
        add_edge!(G, 4i + 2, 4i + 4)
        add_edge!(G, 4i + 3, 4i + 5)
        add_edge!(G, 4i + 4, 4i + 5)
end
add_edge!(G, 4M + 1, 4M + 2)
add_edge!(G, 4M + 2, 4M + 3)
add_edge!(G, 4M + 2, 4M + 4)
add_edge!(G, 4M + 3, 1)
add_edge!(G, 4M + 4, 1)

# Edge lengths
L = EdgeVector(G, 5.0)
for i in 0:M
        L[5i+1] = 25.0
end

# Metric graph
Γ̃ = MetricGraph(G, L)

# discretize MetricGraph
N = 49
Γ, vmap, emap = subdivide_graph(Γ̃, N)

# Initialize
steps = 2^14
u0 = VertexVector(Γ.g, 0.0)
us = Vector{Vector{Float64}}(undef, steps)
for i in 0:2
        # Set initial bumps at first 3 long segments
        u0[vmap[5i+1]] .= 0.5
end

a = 0.2
nagumo(u) = u * (1 - u) * (u - a)

rd_dynamics!(us, u0, nagumo, 2^-7, Γ, vmap, emap, playback=true)
#%%
