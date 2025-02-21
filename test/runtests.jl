#### Files in `test` folder ##
using MetricGraphs
using Test
using StaticArrays
using Graphs

@testset "testing for a certain set of data" begin

    #NOTE: This is temporary
    using Pkg
    using Graphs
    Pkg.activate(".")
    include("./src/Types.jl")
    include("./src/Utils.jl")

    # Laplacian
    a = 1.0
    b = 0.0
    c = 0.0
    H = DivergenceForm(a, b, c)

    # Discrete graph
    G = DiGraph(6)
    add_edge!(G, 1, 2)
    add_edge!(G, 2, 3)
    add_edge!(G, 2, 4)
    add_edge!(G, 4, 6)
    add_edge!(G, 5, 6)
    add_edge!(G, 6, 1)

    # Edge lengths
    L = EdgeFunction(G, 1.0)

    # Metric graph
    Γ = MetricGraph(G,)

    N = 99
    g, vmap, emap = subdivide_graph(G, N)

    steps = 1000
    us = Vector{SVector{nv(g),Float64}}(undef, steps)
    ts = Vector{Float64}(undef, steps)

    # Write your own tests here.
    A = [2.67485 0.186124 0.583946 -0.145416 0.401269;
        -0.17813 0.526082 -0.325439 -1.46627 -3.09161;
        -1.52653 -0.631375 -0.618636 -2.36654 0.291429;
        -0.609349 -1.0533 -0.323403 1.34789 -0.409167;
        0.0687345 -0.946264 -0.118661 -1.49908 -1.30312]

    A = A * A'

    A = (A + A') / 2

    b = [0.38052039048801956;
        0.5554436732038299;
        1.7334160203093987;
        0.4470254916964832;
        1.7676230141050555]

    c = 1.0

    v = ones(5)

    problem_data = testPackageProblem(A, b, c)

    problem_setting = testPackageSetting(γ=2.0)

    result_approx = [0.10055689253917299, 0.24322623669390434, 0.0693976149433379, 0.09597235131516035, -0.530323239453701]

    @test solver_prox_quad(v, problem_data, problem_setting) ≈ result_approx

end

# we can potentially create more matrices like that
