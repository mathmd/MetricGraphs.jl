#%%

using MetricGraphs
using Graphs

#%% Combinatorial Graph G

k = 4 # number of edges on each side of the splitted star graph
g = DiGraph(2k + 2)
add_edge!(g, k + 1, k + 2)
bridge = Edge(k + 1 => k + 2)
for j in 1:k
        add_edge!(g, j, k + 1)
        add_edge!(g, k + 2, 2 + k + j)
end

#%% Visualize metric graph

using Plots, GraphRecipes
gplot = graphplot(g, names=1:nv(g), nodesize=0.3)

#%% Edge Length

l = EdgeVector(g, 25.0) # default length is the truncation of the infinite edges
ℓ = 25.0 # The length of the bridge
l[bridge] = ℓ

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

#%% Initialize

us = Vector{Vector{Float64}}(undef, steps) # empty vector for dynamic solution
u0 = VertexVector(param.Γ.g, 0.0)
for j in 1:k
        # u0[vmap[Edge(j => k+1)][1:N÷3]] .= 1.0 # perturb the first 1/3 on the left edges
        u0[vmap[Edge(j => k + 1)]] .= 1.0 # perturb the whole left edges
end
u0[vmap[bridge][1:2N÷3]] .= 1.0 # purturb the first two half of the bridge

# u0[vmap[bridge][N÷3:2N÷3]] .= 1.0 # perturb middle part of the bridge

#%% Dynamic simulation for symmetric wave from the left

steps = 2^15
rd_dynamics!(us, u0, nagumo, 2^-7, param.Γ, vmap, emap, playback=true)

#%% Define plotting function visualizing solution

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


#%% Define time independent solver and jacobian for Numerical continuation

using BifurcationKit
using LinearAlgebra

function nagumo_ss!(f, x, p, t=0)
        for e in emap[bridge]
                p.Γ.l[e] = p.ℓ / (N + 1)
        end
        L⁻² = Diagonal(1 ./ p.Γ.l .^ 2)

        f .= -Δᵀ * L⁻² * Δᵀ' * x + nagumo(x, p.a)
end

function nagumo_jac(x, p)
        for e in emap[bridge]
                p.Γ.l[e] = p.ℓ / (N + 1)
        end
        L⁻² = Diagonal(1 ./ p.Γ.l .^ 2)

        -Δᵀ * L⁻² * Δᵀ' + Diagonal(x .* (2(p.a + 1) .- 3x) .- p.a)
end

nagumo_ss(x, p, t=0) = nagumo_ss!(similar(x), x, p, t)

#%% 

prob = BifurcationKit.BifurcationProblem(nagumo_ss, us[end], param, (@optic _.ℓ);
        J=nagumo_jac,
        record_from_solution=(x, p; k...) -> (m = x[vmap[bridge][end÷2]]),
        plot_solution=(x, p; kwargs...) -> plot!(x; ylabel="solution", label="", kwargs...),
)

optnewton = NewtonPar(tol=1e-11, verbose=true)
sol = @time BifurcationKit.solve(prob, Newton(), optnewton)

plot_sol(sol.u, param)

#%%

optcont = ContinuationPar(dsmin=0.00001, dsmax=0.2, ds=0.1, p_min=0.0, p_max=25.0,
        newton_options=NewtonPar(max_iterations=10, tol=1e-9), max_steps=1000, n_inversion=4)

br = continuation(prob, PALC(), optcont, bothside=true, normC=norminf; plot=true)

#%% 

plot_sol(br.specialpoint[1].x, param)

#%% bifurcation of 


