	abstract LadderOpFlag
	abstract RaiseOpFlag <: LadderOpFlag
	abstract LowerOpFlag <: LadderOpFlag	

	ladder_inds(n, ::Type{RaiseOpFlag}) = ([1:n-1], [2:n])
	ladder_inds(n, ::Type{LowerOpFlag}) = ([2:n], [1:n-1])
	laddercoeffs(n) = complex(sqrt(linspace(1, n, n)))
	
	function fill_op_arr!(arr::AbstractMatrix, ladderflag)
		if size(arr, 1) == size(arr, 2)
			len = size(arr, 1)
			inds = ladder_inds(len, ladderflag)
			coeffs = laddercoeffs(len)
			for i=1:len-1
				arr[inds[1][i], inds[2][i]] = coeffs[i] 
			end
		else
			error("Cannot generate ladder coefficients for non-square matrix")
		end
		return arr
	end

	# returns a coefficient matrix
	# for a ladder operator for a
	# single particle fock basis
	gen_op_mat(len, ladderflag) = fill_op_arr!(spzeros(Complex128, len, len), ladderflag)

	# this could/should be further optimized,
	# it uses the naive approach of taking the 
	# kronecker product of identity matrices 
	# and the relevant ladder operator matrix
	function laddermatrix(lens, ladderflag, n=1)   
		if n==1
			arr = gen_op_mat(lens[1], ladderflag)
		else
			arr = speye(Complex128, lens[1])
			for i=2:n-1
				arr = kron(speye(Complex128, lens[i]), arr)
			end
			arr = kron(gen_op_mat(lens[n], ladderflag), arr)
		end	
		for i=n+1:length(lens)
			arr = kron(speye(Complex128, lens[i]), arr)
		end
		return arr
	end

	raisematrix(lens, n=1) = laddermatrix(lens, RaiseOpFlag, n)
	lowermatrix(lens, n=1) = laddermatrix(lens, LowerOpFlag, n)

	# n specifies which particle (a.k.a tensor product 
	# factor) this operator acts on
	raiseop(b::AbstractBasis, n=1) = quarr(raisematrix(size(b), n), b, b)
	raiseop(lens::Tuple, n=1) = raiseop(FiniteBasis(lens), n)
	
	lowerop(b::AbstractBasis, n=1) = quarr(lowermatrix(size(b), n), b, b)
	lowerop(lens::Tuple, n=1) = lowerop(FiniteBasis(lens), n)

export raisematrix,
	lowermatrix,
	raiseop,
	lowerop