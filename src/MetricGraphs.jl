module MetricGraphs

using Graphs

include("./Types.jl")
export Hamiltonian, DivergenceForm, MetricGraph

include("./Utils.jl")
export subdivide_graph

using LinearAlgebra

function rd_dynamics!(
    us,
    ts,
    Γ::MetricGraph;
    l₁=1,
    l₂0=1,
    l₃0=1,
    l₄0=1,
    u0=nothing,
    t0=0,
    # Parameters
    # a=0,
    # D=1,
    # # ρ = 0.001 # This is interesting; steady state?
    # # γ = 0.003
    # ρ=0.03,
    # γ=0.1,
    # κ=0,
    source=nothing,
    retrate=nothing,
    nmntr=2^3, # number of snapshots taken during monitoring
    frames=2^7,
    playback=true,
    start=1,
    dt=2^-5,
    adaptive=false,
    label="",
    stop=nothing,
    l₁s,
)

    # define operator for updating the length
    function expL(
        l, # ls = (l₂,l₃,l₄)
        δt,
        dldt, # u = (u|₂,u|₃,u|₄)
    )
        # return max.(1/N,l .+ (δt .* dldt))
        return l .+ (δt .* dldt)
    end
    # elong_rate(u) = (γ .* u)
    function elong_rate(u, γ, β, ρ, θ; l=1, dt=0.01, tol=1 / N)
        dldt = [d for d ∈ (γ .* u .+ β .* u .^ θ) .- ρ]
        dldt .*= (l .+ dt .* dldt .>= tol)
        return (dldt...,)
    end

    # define operator for updating the concentration
    Δᵀ = incidence_matrix(G)
    A = abs.(Δᵀ') ./ 2
    E = max.(0, Δᵀ') # choosing edges for sparse vector of destination vertices
    if a isa Vector
        adv = Diagonal(E * a)
    else
        adv = a .* Diagonal(E * ones(nv(G)))
    end
    λ = ones(nv(G)) # we will store 1/l in this
    # Iₑ = Diagonal(ones(ne(G)))
    Iᵥ = Diagonal(ones(nv(G)))
    # W = sparse(zeros(3, nv(G))) # extrapolation operator
    W = sparse(zeros(2 + num_passive_branches, nv(G))) # extrapolation operator
    # W[1, 1:3] .= [4, -3, 1]
    # W[1, 4N-1:4N+1] .= [1, -3, 4]
    # W[2, 2N-1:2N+1] .= [1, -3, 4]
    # W[3, 3N-1:3N+1] .= [1, -3, 4]
    W[1, 2N-1:2N+1] .= [1, -3, 4]
    W[2, 3N-1:3N+1] .= [1, -3, 4]
    # Δₚ = sparse(zeros(nv(G), 3)) # "difference" operator only for phantom eges
    Δₚ = sparse(zeros(nv(G), 2 + num_passive_branches)) # "difference" operator only for phantom eges
    # Δₚ[1, 1] = -1
    # Δₚ[4N+1, 1] = 1
    # Δₚ[2N+1, 2] = 1
    # Δₚ[3N+1, 3] = 1
    Δₚ[2N+1, 1] = 1
    Δₚ[3N+1, 2] = 1
    ξ = Vector{Float64}(undef, nv(G))
    # ξ[1:N+1] .= 0
    ξ[N+2:2N+1] .= collect(LinRange(1 / N, 1, N))
    ξ[2N+2:3N+1] .= collect(LinRange(1 / N, 1, N))
    # ξ[3N+2:4N+1] .= collect(LinRange(1 / N, 1, N))
    for k ∈ 1:num_passive_branches
        W[2+k, (3+k)N-1:(3+k)N+1] .= [1, -3, 4]
        Δₚ[(3+k)N+1, (2+k)] = 1
        ξ[(2+k)N+2:(3+k)N+1] .= collect(LinRange(1 / N, 1, N))
    end
    Ξ = Diagonal(ξ) # this is non-zero on the branch vertices
    W .*= 1 / 2 # phantom edges at the moving boundaries

    function expT!(
        u, # u∈us
        t,
        δt,
        l,
        dldt
    )
        λ[1:N+1] .= 1 ./ l₁
        λ[N+2:2N+1] .= 1 / l[1] #l[1] is l₂
        λ[2N+2:3N+1] .= 1 / l[2] #l[2] is l₃
        # λ[3N+2:end] .= 1 / l[3]
        ∂ₜl = sparse(zeros(nv(G)))
        ∂ₜl[2N+1] = dldt[1]
        ∂ₜl[3N+1] = dldt[2]
        # ∂ₜl[4N+1] = dldt[3]
        DLDT = sparse(zeros(nv(G)))
        DLDT[N+2:2N+1] .= dldt[1]
        DLDT[2N+2:3N+1] .= dldt[2]
        # DLDT[3N+2:4N+1] .= dldt[3]
        κdldt = zeros(nv(G))
        κdldt[2N+1] = N * κ * dldt[1] / l[1]
        κdldt[3N+1] = N * κ * dldt[2] / l[2]
        # κdldt[4N+1] = N * κ * dldt[3] / l[3]
        for k ∈ 3:2+num_passive_branches
            λ[k*N+2:(k+1)N+1] .= 1 / l[3] #l[3] is l₄
            ∂ₜl[(k+1)N+1] = dldt[3]
            DLDT[k*N+2:(k+1)N+1] .= dldt[3]
            κdldt[(k+1)N+1] = N * κ * dldt[3] / l[3]
        end
        L_inv = Diagonal(E * λ) # living on the edges
        ∂ₜl = Diagonal(∂ₜl)
        DLDT = Diagonal(DLDT)

        β = δt * N .* Diagonal(λ) * (∂ₜl + (-Δᵀ) * (adv * A .- N * D .* L_inv * Δᵀ') - Ξ * DLDT * ((-Δᵀ) * A + Δₚ * W))

        # u .= solve(LinearProblem(Iᵥ .+ β, u .+ (δt .* (s .- κdldt))))
        u .= solve(LinearProblem(Iᵥ .+ β, u .+ (δt .* (source.(ξ, t) .- κdldt))))
        # u .= solve(LinearProblem(Iᵥ .+ β ./ 2, ( Iᵥ - β ./ 2 ) * u .+ (δt .* (s .- κdldt))))

    end

    steps = length(us)

    if u0 === nothing
        # u0 = fill(1 / (l₁ + l₂0 + l₃0 + l₄0), nv(G))
        u0 = fill(1 / (l₁ + l₂0 + l₃0 + num_passive_branches * l₄0), nv(G))
    elseif !isa(u0, Vector)
        u0 = fill(u0, nv(G))
    end

    # if source isa Vector
    #         s = source
    # else
    #         s = zeros(nv(G))
    #         s[1] = source # source at the soma
    # end

    # iteration
    if playback
        nvis = Int(steps / frames)
    else
        nvis = Int(steps / nmntr)
    end

    l2, l3, l4 = l₂0, l₃0, l₄0
    u = copy(u0)
    t = t0
    dldt = 0
    n = 1

    function adjust_dt(l2, l3, l4)
        if adaptive
            return min(1, min(l2, l3, l4)^2)
        else
            return 1
        end
    end

    if isnothing(stop)
        stop = steps
    end

    # leng = sum(l₁[2:end]) / N

    for n ∈ start:stop
        # while l2>0 && l3>0

        δt = dt * adjust_dt(l2, l3, l4)
        t += δt
        # dldt = elong_rate((u[2N+1], u[3N+1], u[4N+1]), γ, β, ρ, θ; l=(l2, l3, l4), dt=δt, tol=1 / N)
        dldt = elong_rate((u[2N+1], u[3N+1], u[4N+1]), γ, β, retrate(t), θ; l=(l2, l3, l4), dt=δt, tol=1 / N)
        # if dl[2]*δt + l₃s[n] < 1/N
        #         dldt = (dl[1],0,dl[3])
        # else
        #         dldt = (dl[1],dl[2],dl[3])
        # end
        l2, l3, l4 = expL((l2, l3, l4), δt, dldt)
        l₂s[n], l₃s[n], l₄s[n] = l2, l3, l4
        expT!(u, t, δt, (l2, l3, l4), dldt)
        # u ./= sum(A * u ./ (E * (1 ./λ))) ./ N
        us[n] = copy(u)
        ts[n] = t
        if n % nvis == 0
            if playback
                # p.edge_color[] = A*u
                # display(heatmap(xc, vc, ρ'; xlims=(xc[1], xc[end]), ylims=(vc[1], vc[end])))
                # plot(
                #         ts[1:n],
                #         l₃s[1:n],
                #         label="l₃"
                # )
                # display(plot!(ts[1:n],l₂s[1:n], label="l₂"))
                plot(ts[1:n], l₁s[1:n] .+ l₂s[1:n], label=label[1])
                plot!(ts[1:n], l₁s[1:n] .+ l₃s[1:n], label=label[3])
                display(plot!(ts[1:n], l₄s[1:n], label=label[2]))
            else
                println("progress: $(100*n/steps) %")
            end
        end
    end
end

export solver_prox_quad

end # module
