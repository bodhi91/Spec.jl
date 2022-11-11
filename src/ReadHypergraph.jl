function ReadHypergraphFile(fname::String)
    i = 0
    j = 0
    num_hyperedges = 0
    num_vertices = 0
    wt_type = 0
    l = 1
    eptr = Int[]
    hedges = Int[]
    hyperedge_wts = Int[]
    vertex_wts = Int[]
    f = open(fname, "r")

    for ln in eachline(f)
        if i == 0
            line_vector = split(ln)
            num_vertices = parse(Int, line_vector[2])
            num_hyperedges = parse(Int, line_vector[1])
            eptr = zeros(Int, num_hyperedges+1)
            eptr[1] = 1
            hyperedge_wts = zeros(Int, num_hyperedges)
            vertex_wts = ones(Int, num_vertices)

            if length(line_vector) > 2
                wt_type = parse(Int, line_vector[3])
            end
        elseif i <= num_hyperedges
            line_vector = split(ln)
            edge_vector = parse.(Int, line_vector)

            if wt_type == 10 || wt_type == 0
                hyperedge_wts[i] = 1
                append!(hedges, edge_vector)
                l += length(edge_vector)
                eptr[i+1] = l
            else
                hyperedge_wts[i] = edge_vector[1]
                append!(hedges, edge_vector[2:end])
                l += length(edge_vector)-1
                eptr[i+1] = l
            end
        else
            vertex_wts[i-num_hyperedges] = parse(Int, ln)
        end
        i += 1
    end

    close(f)

    return (hedges, eptr, vertex_wts, hyperedge_wts, num_vertices, num_hyperedges, length(hedges))
end

function ReadHypergraphFixedFile(fname::String, fixed_vertices::Vector{Int})
    i = 0
    f = open(fname, "r")

    for ln in eachline(f)
        i += 1
        fixed_vertices[i] = parse(Int, ln)
    end
end