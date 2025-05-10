#%%
using MetricGraphs
using Graphs

# Specify the metric graph domain
M = 5
g = DiGraph(4M + 5)
l =
        for i in 0:M-1
                add_edge!(g, 4i + 1, 4i + 2)
                add_edge!(g, 4i + 2, 4i + 3)
                add_edge!(g, 4i + 2, 4i + 4)
                add_edge!(g, 4i + 3, 4i + 5)
                add_edge!(g, 4i + 4, 4i + 5)
        end
add_edge!(g, 4M + 1, 4M + 2)
add_edge!(g, 4M + 2, 4M + 3)
add_edge!(g, 4M + 2, 4M + 4)
# add_edge!(G, 4M + 3, 1)
# add_edge!(G, 4M + 4, 1)
add_edge!(g, 4M + 3, 4M + 5)
add_edge!(g, 4M + 4, 4M + 5)

# Edge lengths
l = EdgeVector(g, 5.0)
ℓ = 18.0
for i in 0:M
        # L[5i+1] = 25.0
        l[5i+1] = ℓ
end

#%% Visualize metric graph

using Plots, GraphRecipes
gplot = graphplot(g, names=1:nv(g), nodesize=0.3)

#%%

# Metric graph
Γ = MetricGraph(g, l)

# discretize MetricGraph
N = 49
# Γ, vmap, emap = subdivide_graph(Γ̃, N)
mgd = MetricGraphDomain(Γ, N)
# Γ, vmap, emap = mgd.Γ̃, mgd.vmap, mgd.emap
# Δᵀ = incidence_matrix(Γ.g)
vertex_embedding = Vector{Tuple{Float64,Float64}}(undef, nv(g))
# for v ∈ 1:num_rungs+2
#         vertex_embedding[v] = (v, 1)
#         vertex_embedding[v+num_rungs+2] = (v, 0)
# end
for m in 0:M
        vertex_embedding[4m+1] = (m, 0)
        vertex_embedding[4m+2] = (m + 0.5, 0.5)
        vertex_embedding[4m+3] = (m + 0.5, -0.5)
        vertex_embedding[4m+4] = (m + 1, 0.0)
end
vertex_embedding[4M+5] = (M + 2, 0.0)
# vertex_embedding[4M+1] = (M,0)
# vertex_embedding[4M+2] = (M+0.5,0.5)
# vertex_embedding[4M+3] = (M+0.5,-0.5)
# vertex_embedding[4M+4] = (M+1,0.0)
embed_metricgraph!(mgd, vertex_embedding)
# Simple visualization
plot(mgd, zeros(nv(mgd.Γ̃.g)), linecolor=:black)

#parameters
param = (Γ=Γ, vmap=vmap, emap=emap, ℓ=ℓ, a=0.2)

# Initialize
steps = 2^15
u0 = VertexVector(Γ.g, 0.0)
us = Vector{Vector{Float64}}(undef, steps)
for i in 0:2
        # Set initial bumps at first 3 long segments
        u0[vmap[5i+1]] .= 0.5
end

nagumo(u, a) = @. u * (1 - u) * (u - a)
nagumo(u) = nagumo(u, param.a)

rd_dynamics!(us, u0, nagumo, mgd, playback=true)

#%% defining parameters
Base.@kwdef mutable struct ModelParameters
        mgd::MetricGraphDomain = mgd # metric graph
        reaction::Function = y -> broadcast(x -> 0.0, y)
        update_length::Union{Function,Nothing} = nothing
        a::Float64 = 0.25 # Nagumo parameter
        ℓ::Float64 = 25.0 # rungs length
        L::Float64 = 25.0 # rail segment length
end
param = ModelParameters(a=0.45)
nagumo(u) = @. u * (1 - u) * (u - param.a)
param.reaction = nagumo
# below is for numerical continuation w.r.t. ℓ
ℓ = 18.0
param.ℓ = ℓ # make the scope of ℓ global
ℓ = 25.0
bridges = []
for i in 0:M
        # L[5i+1] = 25.0
        push!(bridges, collect(edges(mgd.Γ.g))[5i+1])
end
param.update_length() = begin
        for e ∈ bridges
                update_length!(mgd, e, ℓ)
        end
end
#%%

using Plots

function plot_sol(x, p)

        y = 0.0
        flattened_edges = []
        for dedges in emap
                flattened_coordinate = [y]
                for e in dedges
                        y += p.Γ.l[e]
                        push!(flattened_coordinate, y)
                end
                push!(flattened_edges, flattened_coordinate)
        end
        plot(flattened_edges[1], x[vmap[1]], label="", ylim=(0.0, 1.1))
        for i in 2:length(emap)-1
                plot!(
                        flattened_edges[i],
                        x[vmap[i]],
                        label=""
                )
        end
        i = length(emap)
        display(plot!(
                flattened_edges[i],
                x[vmap[i]],
                label=""
        ))
end

#%%

using BifurcationKit
using LinearAlgebra

function nagumo_ss!(f, x, p, t=0)
        for i in 0:M
                for e in p.emap[5i+1]
                        p.Γ.l[e] = p.ℓ / (N + 1)
                end
        end
        L⁻² = Diagonal(1 ./ p.Γ.l .^ 2)

        f .= -Δᵀ * L⁻² * Δᵀ' * x + nagumo(x, p.a)
end

function nagumo_jac(x, p)
        for i in 0:M
                for e in p.emap[5i+1]
                        p.Γ.l[e] = p.ℓ / (N + 1)
                end
        end
        L⁻² = Diagonal(1 ./ p.Γ.l .^ 2)

        -Δᵀ * L⁻² * Δᵀ' + Diagonal(x .* (2(p.a + 1) .- 3x) .- p.a)
end

nagumo_ss(x, p, t=0) = nagumo_ss!(similar(x), x, p, t)

#%%


prob = BifurcationKit.BifurcationProblem(nagumo_ss, us[end], param, (@optic _.ℓ);
        J=nagumo_jac,
        record_from_solution=(x, p; k...) -> (n2=norm(x), s=sum(x), s2=x[end÷2], s4=x[end÷4], s5=x[end÷5]),
        # plot_solution=(x, p; kwargs...) -> (plot!(flattened_edges, x; ylabel="solution", label="", kwargs...))
        # plot_solution=(x, p; kwargs...) -> plot_sol(x, param),
        plot_solution=(x, p; kwargs...) -> plot!(x; ylabel="solution", label="", kwargs...),
        # plot_solution=(x, p; kwargs...) -> (plotsol(x; label="", kwargs...)),
        #          record_from_solution = (x, p; k...) -> x[div(nv(Γ.g),2)]
)

#%%

optnewton = NewtonPar(tol=1e-11, verbose=true)
sol = @time BifurcationKit.solve(prob, Newton(), optnewton)

plot_sol(sol.u, param)

#%%

optcont = ContinuationPar(dsmin=1e-7, dsmax=0.2, ds=0.1, p_min=0.0, p_max=50.0,
        newton_options=NewtonPar(max_iterations=10, tol=1e-9), max_steps=1000, n_inversion=4)

# opts_br_eq = ContinuationPar(dsmin=0.01, dsmax=0.1, ds=0.1,
#         p_min=5.0, p_max=25.0,
#         # newton_options = opt_newton, 
#         # max_steps = 1000,
#         # specific options for precise localization of Hopf points
#         # n_inversion = 6
# )

br = continuation(prob, PALC(), optcont, bothside=true, normC=norminf; plot=true)

#%%

plot_sol(br.specialpoint[1].x, param)

#%% 
