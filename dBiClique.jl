#=function dBiClique(X::AbstractArray, pindex::Pindex)
    n = length(X)
    y = zeros(n)
    d1 = ones(length(pindex.p1))
    d2 = ones(length(pindex.p2))
    t1 = Threads.@spawn (sum(d2) .* d1 .* X[pindex.p1] - d1 * (d2' * X[pindex.p2])) 
    t2 = Threads.@spawn (sum(d1) .* d2 .* X[pindex.p2] - d2 * (d1' * X[pindex.p1])) 
    t1 = fetch(t1)
    t2 = fetch(t2)
    y[pindex.p1] = t1
    y[pindex.p2] = t2

    return y
end=#

function dBiClique(X::AbstractArray, hint_vector::Vector{Vector{Int}})
    n = length(X)
    y = zeros(n)
    d1 = ones(length(hint_vector[1]))
    d2 = ones(length(hint_vector[2]))
    t1 = Threads.@spawn (sum(d2) .* d1 .* X[hint_vector[1]] - d1 * (d2' * X[hint_vector[2]])) 
    t2 = Threads.@spawn (sum(d1) .* d2 .* X[hint_vector[2]] - d2 * (d1' * X[hint_vector[1]])) 
    t1 = fetch(t1)
    t2 = fetch(t2)
    y[hint_vector[1]] = t1
    y[hint_vector[2]] = t2

    return y




    n = length(X)
    y = zeros(n) 
    k = length(hint_vector)
    l = [length(x) for x in hint_vector]
    d = [ones(x) for x in l]
    
    for i in 1:k
        di = d[i]
        for j in i+1:k
            dj = d[j]
            y[hint_vector[i]] += sum(dj) .* di .* X[hint_vector[i]] - di * (dj' * X[hint_vector[j]])
            y[hint_vector[j]] += sum(di) .* dj .* X[hint_vector[j]] - dj * (di' * X[hint_vector[i]])
        end
    end

    return y
end