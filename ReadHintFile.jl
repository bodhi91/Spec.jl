function ReadHintFile(hint::String, new_indices::Vector{Int}, n::Int)
    partition = zeros(Int, n)
    i = 0;
    f = open(hint, "r")

    for ln in eachline(f)
        i += 1

        if new_indices[i] == 0
            continue
        end

        p = parse(Int, ln)
        partition[new_indices[i]] = p
    end

    return partition
end