#TODO: Make N depends on the edges
function subdivide_graph(Γ::MetricGraph, N::Int)
        g = Γ.g
        subdivided = DiGraph()
        vertex_map = Dict{Int,Int}()  # Track mapping of original vertices
        edge_map = Dict{Edge,Vector{Int}}()  # Track mapping for original edges

        # Add all original vertices to the subdivided graph
        for v in vertices(g)
                add_vertex!(subdivided, v)
                vertex_map[v] = v
        end

        # Subdivide each edge
        for e in edges(g)
                u, v = e.src, e.dst  # Source and destination of the edge
                last_vertex = u  # Start at the source vertex
                intermediate_vertices = []

                # Add intermediate vertices and edges
                for i in 1:N
                        new_vertex = nv(subdivided) + 1  # Create a new vertex ID
                        push!(intermediate_vertices, new_vertex)  # Keep track for edge_map
                        add_vertex!(subdivided, new_vertex)
                        add_edge!(subdivided, last_vertex, new_vertex)
                        last_vertex = new_vertex
                end

                # Connect the last intermediate vertex (or source) to the original destination
                add_edge!(subdivided, last_vertex, v)

                # Map the original edge to its intermediate vertices
                edge_map[e] = intermediate_vertices
        end

        return subdivided, vertex_map, edge_map
end
