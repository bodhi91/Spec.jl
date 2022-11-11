function LinearAlgebra.ldiv!(c::AbstractVecOrMat{T}, P::CombinatorialMultigrid.lPreconditioner, b::AbstractVecOrMat{T}) where {T}
	(n, m) = size(b)

    for k in 1:m
        c[:, k] = P.f(b[:, k])
	end
end

@inline function MakeBFunc(vWts::Vector, hint_vector::Vector{Vector{Int}}, multiplier::Vector)
	bfunc = X -> dCliques(X, vWts, multiplier) + 10*dBiClique(X, hint_vector)
	return bfunc
end

#=@inline function MakeBFunc(vWts::Vector, pindex::Pindex, multiplier::Vector)
	bfunc = X -> dCliques(X, vWts, multiplier) + 10*dBiClique(X, pindex)
	return bfunc
end=#

@inline function MakeAFunc(H::Hypergraph)
	afunc = X -> hypL(H, X)
	return afunc
end

function ProjectionStep!(eigvecs::Array{Float64}, degsG::Vector{Float64}, Bfun::Function)
    (n, m) = size(eigvecs)
	o = ones(Float64, n)

	for i in 1:m
		c = -dot(degsG, eigvecs[:, i]) / dot(degsG, o)
		eigvecs[:, i] += c * o
		eigvecs[:, i] /= sqrt(eigvecs[:, i]' * Bfun(eigvecs[:, i]))
	end

	for i in 1:n
		eigvecs[i, :] /= norm(eigvecs[i, :], 2)
	end
end

function GenEigenVecs(H::Hypergraph, A::SparseMatrixCSC, hint_vector::Vector{Vector{Int}}, largest::Bool, nev::Int, solver_iters::Int; bmap = nothing) 
	vWts = H.vwts
	#d = rand(length(vWts)) ./ 1e06
    L = lap(A) #+ spdiagm(d)
	(~, n) = size(L)
	multiplier = ones(size(vWts, 2))
	afunc = MakeAFunc(H)
	amap = LinearMap(afunc, issymmetric=true, n)
	evec = Float64[]

    if bmap == nothing
		bfunc = MakeBFunc(vWts, hint_vector, multiplier)
		bmap = LinearMap(bfunc, n)
	end

    if size(L)[1] < 100
        results = lobpcg(L, largest, nev + 1, maxiter=solver_iters)
        evec = results.X[:,2:nev+1]
    else
        (pfunc, hierarchy) = CombinatorialMultigrid.cmg_preconditioner_lap(L)
        results = lobpcg(amap, bmap, largest, nev, tol=1e-40, maxiter=solver_iters, P=CombinatorialMultigrid.lPreconditioner(pfunc), log=true)
        @info "[EIGEN VECTOR DETAILS] :: $results"
        evec = results.X
    end

	deg = GraphLaplacians.degree_matrix(L)
	deg = deg.diag
	ProjectionStep!(evec, deg, bfunc)

    return evec
end

#=function GenEigenVecs(H::Hypergraph, A::SparseMatrixCSC, pindex::Pindex, largest::Bool, nev::Int, solver_iters::Int; bmap = nothing) 
    vWts = H.vwts
	d = rand(length(vWts)) ./ 1e06
    L = lap(A) + spdiagm(d)
	(~, n) = size(L)
	multiplier = ones(size(vWts, 2))
	Y = rand(n)
	Y = zeros(n)
	Y[pindex.p1] .= -1.0
	Y[pindex.p2] .= 1.0
	afunc = MakeAFunc(H)
	amap = LinearMap(afunc, issymmetric=true, n)
	evec = Float64[]

    if bmap == nothing
		bfunc = MakeBFunc(vWts, pindex, multiplier)
		bmap = LinearMap(bfunc, n)
	end

    if size(L)[1] < 100
        results = lobpcg(L, largest, nev + 1, maxiter=solver_iters)
        evec = results.X[:,2:nev+1]
    else
        (pfunc, hierarchy) = CombinatorialMultigrid.cmg_preconditioner_lap(L)
        results = lobpcg(amap, bmap, largest, nev, tol=1e-40, maxiter=solver_iters, P=CombinatorialMultigrid.lPreconditioner(pfunc), log=true)
        @info "[EIGEN VECTOR DETAILS] :: $results"
        evec = results.X
    end

    return evec
end=#