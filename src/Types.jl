abstract type Hamiltonian end

#TODO: Make a,b,c non-constant
mutable struct DivergenceForm{T<:Union{Real,Array{Real}}} <: Hamiltonian
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

# ✅ Multiple edge keys
# function Base.setindex!(ev::EdgeVector, vals::AbstractVector, es::AbstractVector{<:Union{Edge,Int}})
#         @assert length(vals) == length(es) "Mismatched assignment length"
#         for (e, v) in zip(es, vals)
#                 ev[e] = v
#         end
# end

# --- Vectorized get/set by Edge vector ---
Base.getindex(ev::EdgeVector, es::AbstractVector{<:Edge}) = [ev[e] for e in es]
Base.setindex!(ev::EdgeVector, vals::AbstractVector, es::AbstractVector{<:Edge}) = foreach((e, v) -> ev[e] = v, es, vals)

# Broadcasting into edge index slices
function Base.setindex!(ev::EdgeVector{T}, val::T, es::AbstractVector{<:Edge}) where {T}
        for e in es
                ev[e] = val
        end
end

# Tell Julia that EdgeVector is mutable and supports broadcast setindex!
Base.broadcastable(ev::EdgeVector) = ev

struct MetricGraph{T<:Real}
        g::DiGraph
        l::EdgeVector{T}
        function MetricGraph(g::DiGraph, l::EdgeVector{T}) where {T}
                if length(l) != ne(g)
                        error("Need to specify $(ne(g)) edge lengths.")
                end
                new{T}(g, l)
        end
end

struct MetricGraphDomain{T<:Real}
        Γ::MetricGraph{T}
        res::Int     # number of vertices inserted 
        Γ̃::MetricGraph{T}   # Subdivided graph
        vmap::EdgeVector{Vector{Int}}
        emap::EdgeVector{Vector{Edge}}
        edge_embedding::EdgeVector{Vector{Any}}
        # constructing function
        function MetricGraphDomain(Γ::MetricGraph{T}, res::Int;
                edge_embedding=EdgeVector(Γ.g, [])) where {T}

                Γ̃, vmap, emap = subdivide_graph(Γ, res)
                new{T}(Γ, res, Γ̃, vmap, emap, edge_embedding)
        end
end

@recipe function f(mgd::MetricGraphDomain{T}, u::Vector{T}) where {T<:Real}
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
        zlims --> (0.0, 1.01)
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
