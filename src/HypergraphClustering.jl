function HypergraphCC(H::Hypergraph, excludedEdges::Vector{Int}, mode::Bool)
    eptr = H.eptr
    hedges = H.hedges
    e = H.e
    n = H.n
    id = Vector{Int}(1:n)
    sz = ones(Int, n)

    includedEdges = setdiff(Vector{Int}(1:e), excludedEdges)

    for k in 1:length(includedEdges)
        kk = includedEdges[k]
        l_edge = eptr[kk+1]-eptr[kk]

        for j in 0:l_edge-2
            u = hedges[eptr[kk]+j]
            v = hedges[eptr[kk]+j+1]

            while u != id[u]
                id[u] = id[id[u]]
                u = id[u]
            end

            while v != id[v]
                id[v] = id[id[v]]
                v = id[v]
            end

            if u != v
                if sz[u] < sz[v]
                    id[u] = v
                    sz[v] += sz[u]
                else
                    id[v] = u
                    sz[u] += sz[v]
                end
            end
        end
    end

    cID = zeros(Int, n)
    cs = zeros(Int, n)

    for k in 1:n
        u = k
        
        while u != id[u]
            u = id[u]
        end

        cID[k] = u
        cs[u] += 1
    end

    ndx = findall(x-> x > 0, cs)
    ord = zeros(Int, maximum(ndx))
    ord[ndx] = Vector{Int}(1:length(ndx))

    for k in 1:n
        cID[k] = ord[cID[k]]
    end

    cs = cs[ndx]

    return (cID, cs)
end

function StandardizeClustering(cID::Vector{Int}, cs::Vector{Int})
	cID .+= 2
    cID[1] = 1
    cID[2] = 2
    B = zeros(Int, length(cID))

    @inbounds for j in 1:length(cID)
        B[cID[j]] = 1
    end

    X = findall(!iszero, B)
    M = zeros(Int, maximum(X))
    M[X] = 1:length(X)
    cID = M[cID]
    cs[3:end] = cs[1:end-2]
    cs[1:2] .= 1

    return (cID, cs)
end
