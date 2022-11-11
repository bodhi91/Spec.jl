module Spec

using LDLFactorizations
using Statistics
using IterativeSolvers
using Laplacians
using LinearAlgebra
using LinearMaps
using ParallelKMeans
using Random
using SparseArrays
using GraphLaplacians

# Defining the basic hypergraph structures

mutable struct Hypergraph
    n::Int
    e::Int 
    hedges::Vector{Int}
    eptr::Vector{Int}
    vwts::Vector{Int}
    hwts::Vector{Int}
end

mutable struct HypergraphToken
    hgraph::Hypergraph
    incidence_list::Vector{Vector{Int}}
    hyperedges_pair_list::Array{Int}
    community::Vector{Int}
end

mutable struct Incidence
	n::Int
	e::Int
	hedges::Vector{Int}
	eptr::Vector{Int}
	d::Vector{Int}
end

struct Pindex
    p1::Vector{Int}
    p2::Vector{Int}
end

include("ReadHypergraph.jl")
include("ReadHintFile.jl")
include("cmg/src/CMGAlg.jl")
include("cmg/src/CombinatorialMultigrid.jl")
include("IsolateIslands.jl")
include("HypergraphClustering.jl")
include("HypergraphToGraph.jl")
include("hypL.jl")
include("dCliques.jl")
include("dBiClique.jl")
include("GenerateEmbedding.jl")

export GlobalClusters

Base.@ccallable function GlobalClusters(hg::Ptr{UInt8}, hint::Ptr{UInt8}, Nparts::Cint, nev::Cint, solver_iters::Cint, seed::Cint, expander_cycles::Cint, imbalance::Float32)::Cvoid
    BLAS.set_num_threads(Threads.nthreads())
    
    # Explicit conversion of C types to Julia types
    
    hg = Base.unsafe_string(hg)
    hint = Base.unsafe_string(hint)
    Nparts = Int(Nparts)
    nev = Int(nev)
    solver_iters = Int(solver_iters)
    seed = Int(seed)
    expander_cycles = Int(expander_cycles)
    imbalance = Float64(imbalance)
    
    println("Starting Spectral Interface")

    hg_split = split(hg, "/")
    hg_name = hg_split[end]
    hypergraph_file = hg
    num_parts = Nparts
    ub_factor = imbalance
    Random.seed!(seed)

    line_log = repeat("=", 80)
    @info "$line_log"
    println("STARTING SUPERVISED SPECTRAL PARTITIONING ENGINE")
    @info "$line_log"

    (hedges, eptr, vertex_weights, hyperedge_weights, num_vertices, num_hyperedges, ~) = ReadHypergraphFile(hypergraph_file)
    hypergraph = Hypergraph(num_vertices, num_hyperedges, hedges, eptr, vertex_weights, hyperedge_weights)
    hsizes = hypergraph.eptr[2:end] - hypergraph.eptr[1:end-1]
    hsize_max = maximum(hsizes)
    (processed_hypergraph, original_indices, new_indices, unused_indices) = IsolateIslands(hypergraph)
    
    max_size = sum(vertex_weights) * (50+ub_factor)/100
    total_vwts = sum(vertex_weights)
    hint_partition = ReadHintFile(hint, new_indices, processed_hypergraph.n)
    #hint_partition = -ones(Int, processed_hypergraph.n)
    
    println("TOTAL VERTICES: $num_vertices")
    println("TOTAL HYPEREDGES: $num_hyperedges")
    println("POST PROCESSING :: TOTAL VERTICES: $(processed_hypergraph.n)")
    println("POST PROCESSING :: TOTAL HYPEREDGES: $(processed_hypergraph.e)")
    println("SIZE OF LARGEST HYPEREDGE: $hsize_max")
    println("MAX CAPACITY CONSTRAINT: $max_size")

    adj_mat = HypergraphToGraph(processed_hypergraph, expander_cycles)

    hint_vector = [Vector{Int}() for i in 1:num_parts]

    for i in 1:length(hint_partition)
        push!(hint_vector[hint_partition[i]+1], i)
    end

    eigenvecs = GenEigenVecs(processed_hypergraph, adj_mat, hint_vector, false, nev, solver_iters)
    eigenvecs_matrix = eigenvecs'
    clusters_token = kmeans(Hamerly(), eigenvecs_matrix, num_parts)
    remapped_eigenvecs = zeros(hypergraph.n, nev)
    clusters_assignments = clusters_token.assignments
    remapped_clusters_assignments = zeros(Int, hypergraph.n)
    centroid = mean(eigenvecs, dims = 1)
    centroid = centroid[1, :]

    f = open("eigenvecs.dat", "w")

    for i in 1:hypergraph.n
	for j in 1:nev
	    print(f, remapped_eigenvecs[i, j], " ")
	end

	println(f, "\n")
    end

    close(f)

    new_idx::Int = 0

    for i in 1:length(remapped_clusters_assignments)
        new_idx = new_indices[i]
        
        if new_idx == 0 
            continue
        end
        
	remapped_eigenvecs[i, :] = eigenvecs[new_idx, :]
        remapped_clusters_assignments[i] = clusters_assignments[new_idx]
    end

    for i in 1:length(unused_indices)
	remapped_eigenvecs[i, :] = centroid;
        remapped_clusters_assignments[unused_indices[i]] = rand(1:num_parts)
    end

    println("Finished Generated K means")
    
    f = open("clusters.dat", "w")

    for x in clusters_assignments
        println(f, x)
    end

    close(f)

    f = open("eigenvecs.dat", "w")

    for i in 1:hypergraph.n
	for j in 1:nev
	    print(f, remapped_eigenvecs[i, j], " ")
	end

	print(f, "\n")
    end

    close(f)
end

end #end module
