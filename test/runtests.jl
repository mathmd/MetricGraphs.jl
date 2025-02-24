#### Files in `test` folder ##
using MetricGraphs
using Test
using Graphs
using LinearAlgebra
using LinearSolve
using SparseArrays
using Plots, LaTeXStrings

@testset "Testing" begin

    #NOTE: This is temporary
    using Pkg
    Pkg.activate("./test")
    using Graphs
    using StaticArrays
    include("../src/Types.jl")
    include("../src/Utils.jl")

    # Laplacian
    a = 1.0
    b = 0.0
    H = DivergenceForm(a, b)

    #%% Periodic graph
    M = 4
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
    # add_edge!(G, 4M+3, 1)
    # add_edge!(G, 4M+4, 1)

    # Edge lengths
    L = EdgeVector(G, 1.0)
    for i in 0:M
        L[5i+1] = 25.0
    end

    # δ-type Vertex conditions
    α = VertexVector(G, 1.0) # α=0 means Dirichlet condition
    ζ = VertexVector(G, 0.0) # γ=0 means Neumann-Kirchoff condition

    # Metric graph
    Γ̃ = MetricGraph(G, L)

    N = 49
    Γ, vmap, emap = subdivide_graph(Γ̃, N)

    steps = 2^10
    u0 = VertexVector(Γ.g, 0.0)
    us = Vector{Vector{Float64}}(undef, steps)
    for i in 1:11
        u0[vmap[i]] .= 1.0
    end


    #TODO: make reaction term depends on time and space, and include source in it. 
    nagumo(u, a) = u * (1 - u) * (u - a)
    a = 0.2
    reaction(u) = nagumo(u, a)

    rd_dynamics!(us, u0, Γ, vmap, emap)

    #%%

    function rd_dynamics!(
        us,
        u0,
        Γ::MetricGraph,
        vmap,
        emap;
        nmntr=2^3, # number of snapshots taken during monitoring
        frames=2^7,
        playback=true,
        start=1,
        δt=2^-5,
        stop=nothing,
        ylims=(0.0, 1.0)
    )

        # define operator for updating the concentration
        Δᵀ = incidence_matrix(Γ.g)
        # A = abs.(Δᵀ') ./ 2

        Iᵥ = Diagonal(ones(nv(Γ.g)))

        L_inv = Diagonal(1 ./ Γ.l .^ 2) # living on the edges

        function expT!(
            u # u∈us
        )

            β = δt * ((-Δᵀ) * L_inv * (-L_inv * Δᵀ'))

            u .= solve(LinearProblem(Iᵥ .+ β, u .+ (δt .* (reaction.(u)))))
        end

        steps = length(us)

        # iteration
        if playback
            nvis = Int(steps / frames)
        else
            nvis = Int(steps / nmntr)
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

        if isnothing(stop)
            stop = steps
        end

        for n ∈ start:stop

            expT!(u)
            us[n] = copy(u)

            if n % nvis == 0
                if playback
                    plot(flattened_edges[1], u[vmap[1]], ylims=ylims, label="")
                    for i in 2:length(emap)-1
                        plot!(
                            flattened_edges[i],
                            u[vmap[i]],
                            label=""
                            # raw"$e^"*string(i)*raw"$"
                        )
                    end
                    i = length(emap)
                    display(plot!(
                        flattened_edges[i],
                        u[vmap[i]],
                        label=""
                        # raw"$e^{"*string(i)*raw"}$"
                    )
                    )
                else
                    println("progress: $(100*n/steps) %")
                end
            end
        end
    end



    # @test solver_prox_quad(v, problem_data, problem_setting) ≈ result_approx

end

# we can potentially create more matrices like that
