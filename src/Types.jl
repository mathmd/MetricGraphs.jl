abstract type Hamiltonian end

#TODO: Make a,b,c non-constant
mutable struct DivergenceForm{T<:Union{Float64,Array{Float64}}} <: Hamiltonian
        a::T
        b::T
end

struct EdgeVector{T} <: AbstractVector{T}
        g::DiGraph
        values::Dict{Edge,T}  # Map from Edge => Value
        default::T  # Default value for edges not in Dict

        function EdgeVector(g::DiGraph, default::T) where {T}
                values = Dict{Edge,T}()
                new{T}(g, values, default)
        end
end

# Define Base methods for compatibility
Base.size(f::EdgeVector) = (ne(f.g),)  # Number of edges

# Indexing via Edge object
function Base.getindex(f::EdgeVector, e::Edge)
        get(f.values, e, f.default)  # Return stored value or default
end

# Indexing via integer (i.e., using edges(g))
function Base.getindex(f::EdgeVector, i::Int)
        e = collect(edges(f.g))[i]  # Retrieve edge at index i
        return f[e]
end

# Set value via Edge object
function Base.setindex!(f::EdgeVector, value, e::Edge)
        f.values[e] = value
end

# Set value via integer index
function Base.setindex!(f::EdgeVector, value, i::Int)
        e = collect(edges(f.g))[i]  # Retrieve edge at index i
        f[e] = value
end

struct MetricGraph
        g::DiGraph
        l::EdgeVector{Float64}
        function MetricGraph(g, l)
                if length(l) != ne(g)
                        error("Need to specify $(ne(g)) edge lengths.")
                end
                new(g, l)
        end
end

struct MetricGraphDomain
        Γ::MetricGraph
        res::Int     # number of vertices inserted 
        Γ̃::MetricGraph   # Subdivided graph
        vmap::EdgeVector{Vector{Int}}
        emap::EdgeVector{Vector{Edge}}
        edge_embedding::EdgeVector{Vector{Any}}
        # constructing function
        function MetricGraphDomain(Γ::MetricGraph, res::Int; edge_embedding=EdgeVector(Γ.g, []))
                Γ̃, vmap, emap = subdivide_graph(Γ, res)
                new(Γ, res, Γ̃, vmap, emap, edge_embedding)
        end
end

@recipe function f(mgd::MetricGraphDomain, u::Vector{Float64})
        # Plot the metric graph embedded into the xy plane
        @series begin
                seriestype := :path3d
                linecolor := :black
                label := ""
                [(mgd.edge_embedding[e]..., zeros(length(mgd.vmap[e]))) for e in edges(mgd.Γ.g)]
        end
        # Plot the function u on metric graph
        @series begin
                seriestype := :path3d
                label := ""
                linecolor --> :blue
                [(mgd.edge_embedding[e]..., u[mgd.vmap[e]]) for e in edges(mgd.Γ.g)]
        end
        linecolor --> :blue
        label --> ""
        aspect_ratio --> :equal
        zlims --> (0.0, 1.0)
        []
end

# Vertex-based function vector
#WARNING: Right now, deleting vertices mess up edge ordering; this structure is only good for adding new vertices. 
struct VertexVector{T} <: AbstractVector{T}
        g::DiGraph
        values::Dict{Int,T}
        default::T

        function VertexVector(g::DiGraph, default::T) where {T}
                values = Dict{Int,T}()
                new{T}(g, values, default)
        end
end

Base.size(vv::VertexVector) = (nv(vv.g),)

# Indexing via Vertex
function Base.getindex(vv::VertexVector, v::Int)
        get(vv.values, v, vv.default)
end

# Set value via Vertex
function Base.setindex!(vv::VertexVector, value, v::Int)
        vv.values[v] = value
end
