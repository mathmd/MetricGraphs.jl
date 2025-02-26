module MetricGraphs

using Graphs, LinearAlgebra, LinearSolve, SparseArrays
using Plots
using ProgressBars

include("./Types.jl")
export MetricGraph, VertexVector, EdgeVector

include("./Utils.jl")
export subdivide_graph

#TODO: make reaction term depends on time and space, and include source in it. 
function rd_dynamics!(
        us,
        u0,
        reaction,
        δt,
        Γ::MetricGraph,
        vmap,
        emap;
        frames=2^7,
        playback=false,
        start=1,
        stop=nothing,
        ylims=(0.0, 1.0)
)

        Δᵀ = incidence_matrix(Γ.g)
        # A = abs.(Δᵀ') ./ 2

        Iᵥ = Diagonal(ones(nv(Γ.g)))

        L⁻² = Diagonal(1 ./ Γ.l .^ 2)
        β = Iᵥ .+ δt * (Δᵀ * L⁻² * Δᵀ')

        #TODO: Implement second order accuracy, and optimize solving LinearProblem.
        function expT!(u)
                u .= solve(LinearProblem(β, u .+ (δt .* reaction.(u))))
        end

        x = 0.0
        flattened_edges = []
        for dedges in emap
                flattened_coordinate = [x]
                for e in dedges
                        x += Γ.l[e]
                        push!(flattened_coordinate, x)
                end
                push!(flattened_edges, flattened_coordinate)
        end

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
                        plot(flattened_edges[1], u[vmap[1]], ylims=ylims, label="")
                        for i in 2:length(emap)-1
                                plot!(
                                        flattened_edges[i],
                                        u[vmap[i]],
                                        label=""
                                )
                        end
                        i = length(emap)
                        display(plot!(
                                flattened_edges[i],
                                u[vmap[i]],
                                label=""
                        ))
                end
        end
end

export rd_dynamics!

end # module
