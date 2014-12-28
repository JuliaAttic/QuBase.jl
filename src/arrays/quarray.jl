import Base: size,
	length,
	getindex,
	similar,
	in,
	#TODO: Implement the below operations
	*,.*,
	/,./,
	+,.+,
	-,.-,
	kron,
	ctranspose,
	transpose

abstract AbstractQuArray{B<:(AbstractBasis...), T, N} <: AbstractArray{T,N}

###########
# QuArray # 				
###########

    checksize(::Type{Ket}, qa) = ndims(qa) <= 2 && size(qa, 2) == 1 
    checksize(::Type{Bra}, qa) = ndims(qa) <= 2 && size(qa, 1) == 1 

	# this is done separately from the N-dimensional
	# check to handle the fact that row vectors are
	# N=2
	function checkbases(coeffs, bases::(AbstractBasis,))
		if ndims(coeffs) <= 2
			if size(coeffs, 2) == 1
				return checkcoeffs(coeffs, 1, first(bases))
			elseif size(coeffs, 1) == 1
				return checkcoeffs(coeffs, 2, first(bases))
			end	
		end
		return false
	end

	function checkbases(coeffs, bases::(AbstractBasis...))
		if ndims(coeffs) == length(bases)
			return reduce(&, [checkcoeffs(coeffs, i, bases[i]) for i=1:ndims(coeffs)])
		end
		return false
	end

	type QuArray{B<:(AbstractBasis...), T, N, A} <: AbstractQuArray{B, T, N}
	    coeffs::A
	    bases::B
	    function QuArray(coeffs::AbstractArray{T, N}, bases::B) 
	    	if checkbases(coeffs, bases) 
	    		new(coeffs, bases)
	    	else 
	    		error("Coefficient array does not conform to input bases")
	    	end
	    end
	end

	# We can't use N to tell whether something is really a vector
	# since row vectors are N=2 (quite unfortunate). We can, however,
	# infer the behavior from the number of bases...
	typealias QuVector{B<:AbstractBasis,T,N,A} QuArray{(B,),T,N,A}
	typealias QuMatrix{R<:AbstractBasis,C<:AbstractBasis,T,N,A} QuArray{(R,C),T,N,A}
	typealias QuTensor{B<:(AbstractBasis...),T,N,A} QuArray{B,T,N,A}

	QuArray{T,N}(coeffs::AbstractArray{T, N}, bases::Tuple) = QuArray{typeof(bases), T, N, typeof(coeffs)}(coeffs, bases)
	QuArray(coeffs, bases...) = QuArray(coeffs, bases)

	getbasis(qa::QuArray, i) = qa.bases[i]
	getcoeffs(qa::QuArray) = qa.coeffs

	########################
	# Array-like functions #
	########################
	size(qa::QuArray, i...) = size(qa.coeffs, i...)

	similar{B,T}(qa::QuArray{B,T}, element_type, dims) = quarr(similar(qa.coeffs, T, dims), makebasis(dims))
	similar{B,T}(qa::QuArray{B,T}, element_type=T) = quarr(similar(qa.coeffs, T), qa.bases)

	getindex(qa::QuArray, i::AbstractArray) = getindex(qa.coeffs, i)
	getindex(qa::QuArray, i::Real) = getindex(qa.coeffs, i)
	getindex(qa::QuArray, i) = getindex(qa.coeffs, i)
	getindex(qa::QuArray, i...) = getindex(qa.coeffs, i...)

	in(c, qa::QuArray) = in(c, qa.coeffs)

	ctranspose(qa::QuVector) = QuArray(ctranspose(qa.coeffs), qa.bases)
	ctranspose(qa::QuMatrix) = QuArray(ctranspose(qa.coeffs), reverse(qa.bases))
	transpose(qa::QuVector) = QuArray(transpose(qa.coeffs), qa.bases)
	transpose(qa::QuMatrix) = QuArray(transpose(qa.coeffs), reverse(qa.bases))

	######################
	# Include Statements #
	######################
	include("diracarrays.jl")
	include("ladderops.jl")

	############################
	# Convenience Constructors #
	############################
	makebasis(lens::Tuple, B::DataType=FiniteBasis) = ntuple(length(lens), n->B(lens[n]))

	quarr(coeffs) = QuArray(coeffs, makebasis(size(coeffs), FiniteBasis))
	quarr(coeffs, bases::(AbstractBasis...)) = QuArray(coeffs, bases)
	quarr(coeffs, bases::AbstractBasis...) = quarr(coeffs, bases)

	one_at_ind!(arr, i) = setindex!(arr, one(eltype(arr)), i)
	single_coeff(i, lens...) = one_at_ind!(zeros(Complex128, lens), i)
	diraccoeffs(i, len, ::Type{Ket}) = single_coeff(i, len)
	diraccoeffs(i, len, ::Type{Bra}) = single_coeff(i, 1, len)

	statevec(i::Int, fb::FiniteBasis) = quarr(single_coeff(i, length(fb)), fb)
	statevec(i::Int, lens::Int...=s) = statevec(i, FiniteBasis(lens))

	diracvec(coeffs::AbstractArray, D=Ket, S=AbstractStructure) = DiracVector(coeffs, FockBasis{S}(length(coeffs)), D)
	diracvec(i::Int, b::AbstractLabelBasis, D=Ket) = DiracVector(diraccoeffs(i, length(b), D), b, D)
	diracvec(tup::(Int...), b::AbstractLabelBasis, D=Ket) = DiracVector(diraccoeffs(getpos(b, tup), length(b), D), b, D)

	ketvec(s::(Int...), lens::Int...) = diracvec(s, FockBasis(lens), Ket)
	ketvec(s::(Int...)) = diracvec(s, FockBasis(map(x->x+1, s)), Ket)
	ketvec(s::Int, lens::Int...=s) = diracvec(s, FockBasis(lens), Ket)

	bravec(s::(Int...), lens::Int...) = diracvec(s, FockBasis(lens), Bra)
	bravec(s::(Int...)) = diracvec(s, FockBasis(map(x->x+1, s)), Bra)
	bravec(s::Int, lens::Int...=s) = diracvec(s, FockBasis(lens), Bra)

export AbstractQuArray,
	QuArray,
	quarr,
	statevec,
	ketvec,
	bravec
