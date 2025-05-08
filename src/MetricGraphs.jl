module MetricGraphs

using Graphs, LinearAlgebra, LinearSolve, SparseArrays
using Plots
using ProgressBars

include("./Types.jl")
export MetricGraph, VertexVector, EdgeVector, MetricGraphDomain

include("./Utils.jl")
export subdivide_graph, update_length!, embed_metricgraph!

#TODO: make reaction term depends on time and space, and include source in it. 
function rd_dynamics!(
        us,
        u0,
        # reaction,
        # Γ::MetricGraph,
        # vmap,
        # emap;
        # mgd;
        p;
        δt=2^-7,
        frames=2^7,
        playback=true,
        start=1,
        stop=nothing,
        ylims=(0.0, 1.0)
)
        Γ = p.mgd.Γ̃
        vmap = p.mgd.vmap
        emap = p.mgd.emap

        Δᵀ = incidence_matrix(Γ.g)
        # A = abs.(Δᵀ') ./ 2

        Iᵥ = Diagonal(ones(nv(Γ.g)))

        L⁻² = Diagonal(1 ./ Γ.l .^ 2)
        β = Iᵥ .+ δt * (Δᵀ * L⁻² * Δᵀ')

        #TODO: Implement second order accuracy, and optimize solving LinearProblem.
        function expT!(u)
                u .= solve(LinearProblem(β, u .+ (δt .* p.reaction(u))))
        end

        # x = 0.0
        # flattened_edges = []
        # for dedges in emap
        #         flattened_coordinate = [x]
        #         for e in dedges
        #                 x += Γ.l[e]
        #                 push!(flattened_coordinate, x)
        #         end
        #         push!(flattened_edges, flattened_coordinate)
        # end

        u = Vector{Float64}(u0)
        n = 1

        steps = length(us)

        if isnothing(stop)
                stop = steps
        end

        if playback
                nvis = floor(stop / frames)
        end

        for n ∈ ProgressBar(start:stop)

                expT!(u)
                #TODO: Get rid of copy
                us[n] = copy(u)

                if playback && n % nvis == 0
                        #         plot(flattened_edges[1], u[vmap[1]], ylims=ylims, label="")
                        #         for i in 2:length(emap)-1
                        #                 plot!(
                        #                         flattened_edges[i],
                        #                         u[vmap[i]],
                        #                         label=""
                        #                 )
                        #         end
                        #         i = length(emap)
                        #         display(plot!(
                        #                 flattened_edges[i],
                        #                 u[vmap[i]],
                        #                 label=""
                        #         ))
                        display(plot(p.mgd, u))
                end
        end
end

function graph_laplacian(Γ::MetricGraph)
        Δᵀ = incidence_matrix(Γ.g)
        L⁻² = Diagonal(1 ./ Γ.l .^ 2)

        Δᵀ * L⁻² * Δᵀ'
end

function rd_ss!(f, x, p, t=0)
        p.update_length()
        L⁻² = Diagonal(1 ./ p.mgd.Γ̃.l .^ 2)
        Δᵀ = incidence_matrix(p.mgd.Γ.g)

        f .= -Δᵀ * L⁻² * Δᵀ' * x + p.reaction(x, p.a)
        f .= p.Δ * x + p.f(x)
end

function nagumo_jac(x, p)
        p.update_length()
        L⁻² = Diagonal(1 ./ p.mgd.Γ̃.l .^ 2)
        Δᵀ = incidence_matrix(p.mgd.Γ.g)
        -Δᵀ * L⁻² * Δᵀ' + Diagonal(x .* (2(p.a + 1) .- 3x) .- p.a)
end

rd_ss(x, p, t=0) = rd_ss!(similar(x), x, p, t)

export rd_dynamics!, rd_ss!, rd_ss, graph_laplacian, nagumo_jac

end # module
