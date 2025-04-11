#%%

using MetricGraphs
using Graphs

#%% Combinatorial Graph G

k = 4 # number of edges on each side of the splitted star graph
g = DiGraph(2k + 2)
add_edge!(g, 1, 2)
for j in 1:k
        add_edge!(g, 2 + j, 1)
        add_edge!(g, 2, 2 + k + j)
end

#%% Edge Length

l = EdgeVector(g, 25.0) # default length is the truncation of the infinite edges
ℓ = 10.0 # The length of the bridge
l[Edge(1 => 2)] = ℓ

#%% Metric Graph

Γ = MetricGraph(g, l)

#%% discretization

N = 49
Γ̃, vmap, emap = subdivide_graph(Γ, N)
Δᵀ = incidence_matrix(Γ̃.g)

#%% parameters

param = (Γ=Γ̃, vmap=vmap, emap=emap, ℓ=ℓ, a=0.4)
nagumo(u, a) = @. u * (1 - u) * (u - a)
nagumo(u) = nagumo(u, param.a)

#%% Dynamic simulation for symmetric wave from the left

steps = 2^15
us = Vector{Vector{Float64}}(undef, steps) # empty vector for dynamic solution
u0 = VertexVector(param.Γ.g, 0.0)
for j in 1:k
        u0[vmap[Edge(2 + j => 1)][1:N÷3]] .= 1.0 # perturb the first 1/3 on the left edges
end
rd_dynamics!(us, u0, nagumo, 2^-7, param.Γ, vmap, emap, playback=true)

#%%
