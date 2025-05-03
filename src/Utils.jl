#TODO: Make N depends on the edges
#TODO: Make subdivide_graph into a structure so that it automatically subdivide Γ when the main graph changes.

function update_length!(mgd::MetricGraphDomain, e::Union{Edge,Int}, new_length::Float64)
        mgd.Γ.l[e] = new_length
        for edge ∈ mgd.emap[e]
                mgd.Γ̃.l[edge] = new_length / (mgd.res + 1)
        end
end

function subdivide_graph(Γ::MetricGraph, N::Int)
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
