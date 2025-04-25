#%%
using MetricGraphs
using Graphs

# Specify the metric graph domain
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
# add_edge!(G, 4M + 3, 1)
# add_edge!(G, 4M + 4, 1)
add_edge!(G, 4M + 3, 4M + 5)
add_edge!(G, 4M + 4, 4M + 5)

# Edge lengths
L = EdgeVector(G, 5.0)
ℓ = 18.0
for i in 0:M
        # L[5i+1] = 25.0
        L[5i+1] = ℓ
end

#%% Visualize metric graph

using Plots, GraphRecipes
gplot = graphplot(G, names=1:nv(G), nodesize=0.3)

#%%

# Metric graph
Γ̃ = MetricGraph(G, L)

# discretize MetricGraph
N = 49
Γ, vmap, emap = subdivide_graph(Γ̃, N)
Δᵀ = incidence_matrix(Γ.g)

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

rd_dynamics!(us, u0, nagumo, 2^-7, Γ, vmap, emap, playback=true)

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
