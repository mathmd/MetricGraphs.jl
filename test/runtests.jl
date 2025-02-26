#### Files in `test` folder ##
using MetricGraphs
using Test
using Graphs

@testset "Manufacturing solution" begin

    #%%
    u_exact = (x, t) -> exp.(-t) * sin.(x)

    # circle
    G = DiGraph(3)
    add_edge!(G, 1, 2)
    add_edge!(G, 2, 3)
    add_edge!(G, 3, 1)

    # Edge lengths
    L = EdgeVector(G, 2π / 3)

    # δ-type Vertex conditions
    # α = VertexVector(G, 1.0) # α=0 means Dirichlet condition
    # ζ = VertexVector(G, 0.0) # γ=0 means Neumann-Kirchoff condition

    # Metric graph
    Γ̃ = MetricGraph(G, L)

    N = 999
    M = 3(N + 1)
    Γ, vmap, emap = subdivide_graph(Γ̃, N)

    steps = 100
    δt = 2^-9
    u0 = VertexVector(Γ.g, 0.0)
    us = Vector{Vector{Float64}}(undef, steps)
    u0[unique(reduce(vcat, vmap))] .= sin.(range(0.0, 2π * (1 - 1 / M), length=M))

    reaction = u -> 0.0

    rd_dynamics!(us, u0, reaction, δt, Γ, vmap, emap)

    @test length(unique(reduce(vcat, vmap))) == M

    #TODO: error is still quite big; need to do higher order approximation
    @test isapprox(sum((us[end][unique(reduce(vcat, vmap))] .- u_exact(range(0.0, 2π * (1 - 1 / M), length=M), steps * δt)) .^ 2) / M, 0.0, atol=1e-7)
    #%%

end
