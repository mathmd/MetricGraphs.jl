abstract type Hamiltonian end

#TODO: Make a,b,c non-constant
mutable struct DivergenceForm{T<:Union{Float64,Array{Float64}}} <: Hamiltonian
        a::T
        b::T
        c::T
end

struct EdgeVector{T} <: AbstractVector{T}
        g::SimpleDiGraph
        values::Dict{Edge,T}  # Map from Edge => Value
        default::T  # Default value for edges not in Dict

        function EdgeVector(g::SimpleDiGraph, default::T) where {T}
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

