#TODO: Make N depends on the edges
#TODO: Make subdivide_graph into a structure so that it automatically subdivide Γ when the main graph changes.

function update_length!(mgd::MetricGraphDomain{T}, e::Union{Edge,Int}, new_length::Real) where {T<:Real}
        """
        Update length of an edge `e` in `mgd.Γ` (`mgd.Γ.l[e]`) with `new_length` along with the corresponding discritization `mgd.Γ̃`.
        """
        mgd.Γ.l[e] = new_length
        for edge ∈ mgd.emap[e]
                mgd.Γ̃.l[edge] = new_length / (mgd.res + 1)
        end
end

function update_length!(mgd::MetricGraphDomain{T}, edges::Union{Vector{Edge},Vector{Int}}, new_length::Real) where {T<:Real}
        """
        Update multiple edges in `edges`
        """
        for e ∈ edges
                update_length!(mgd, e, new_length)
        end
end

function update_length!(lengths::EdgeVector{T}, mgd::MetricGraphDomain, e::Union{Edge,Int}, new_length::T) where {T<:Real}
        # lengths = EdgeVector{T}(mgd.Γ̃.g, T(0))
        for edge ∈ mgd.emap[e]
                lengths[edge] = new_length / T(mgd.res + 1)
        end
end

function embed_metricgraph!(mgd::MetricGraphDomain, v_embed::Vector{Tuple{T,T}}) where {T<:Real}
        for e ∈ edges(mgd.Γ.g)
                mgd.edge_embedding[e] = []
                o, t = v_embed[src(e)], v_embed[dst(e)]
                push!(mgd.edge_embedding[e], range(o[1], t[1], mgd.res + 2), range(o[end], t[end], mgd.res + 2))
        end
end

function subdivide_graph(Γ::MetricGraph{T}, N::Int) where {T<:Real}
        G = Γ.g
        L = Γ.l
        subdivided = DiGraph(nv(G))
        # vertex_map = Dict{Edge,Vector{Int}}()  # Track mapping of original vertices
        vertex_map = EdgeVector(G, Int[])  # Track mapping of original vertices
        # edge_map = Dict{Edge,Vector{Edge}}()  # Track mapping for original edges
        edge_map = EdgeVector(G, Edge[])  # Track mapping for original edges
        l = EdgeVector(subdivided, 1.0)

        # Add all original vertices to the subdivided graph
        # for v in vertices(g)
        #         add_vertex!(subdivided)
        #         vertex_map[v] = v
        # end

        # Subdivide each edge
        for e in edges(G)
                u, v = e.src, e.dst  # Source and destination of the edge
                last_vertex = u  # Start at the source vertex
                vertices_on_e = [u]
                edges_on_e = []
                λ = L[e] / (N + 1)

                # Add intermediate vertices and edges
                for i in 1:N
                        new_vertex = nv(subdivided) + 1  # Create a new vertex ID
                        add_vertex!(subdivided)
                        small_e = Edge(last_vertex => new_vertex)
                        add_edge!(subdivided, small_e)
                        l[small_e] = λ
                        push!(vertices_on_e, new_vertex)  # Keep track for edge_map
                        push!(edges_on_e, small_e)  # Keep track for edge_map
                        last_vertex = new_vertex
                end

                # Connect the last intermediate vertex (or source) to the original destination
                small_e = Edge(last_vertex => v)
                add_edge!(subdivided, small_e)
                l[small_e] = λ
                push!(vertices_on_e, v)  # Keep track for edge_map
                push!(edges_on_e, small_e)  # Keep track for edge_map

                # Map the original edge to its intermediate vertices
                edge_map[e] = edges_on_e
                vertex_map[e] = vertices_on_e
        end

        return MetricGraph(subdivided, l), vertex_map, edge_map
end

# function plot(mgd::MetricGraphDomain,u::Vector{Float64},kwargs...)
#         # Plot the metric graph embedded into the xy plane
#         plot(,kwargs...)
#         # Plot u on each edge of mgd.Γ
#         for e ∈ edges(mgd.Γ)
#                 v_on_e = mgd.vmap[e]
#         end
# end
