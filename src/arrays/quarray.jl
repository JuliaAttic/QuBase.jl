import Base: 
	size,
	length,
	getindex

abstract AbstractQuArray{B<:(AbstractBasis...), T, N} <: AbstractArray{T,N}
typealias AbstractQuVector{B<:(AbstractBasis,), T} AbstractQuArray{B, T}
typealias AbstractQuMatrix{B<:(AbstractBasis, AbstractBasis), T} AbstractQuArray{B, T}

###########
# QuArray # 				
###########

	checkcoeffs(coeffs, bases::(AbstractBasis...)) = reduce(&, [checkcoeffs(coeffs, i, bases[i]) for i=1:ndims(coeffs)])

	type QuArray{B<:(AbstractBasis...), T, N, A} <: AbstractQuArray{B, T, N}
	    coeffs::A
	    bases::B
	    function QuArray(coeffs::AbstractArray{T, N}, bases::B) 
	    	if checkcoeffs(coeffs, bases) 
	    		new(coeffs, bases)
	    	else 
	    		error("Coefficient array does not conform to input bases")
	    	end
	    end
	end

	QuArray{T,N}(coeffs::AbstractArray{T, N}, bases::Tuple) = QuArray{typeof(bases), T, N, typeof(coeffs)}(coeffs, bases)
	QuArray(coeffs, bases...) = QuArray(coeffs, bases)

	########################
	# Array-like functions #
	########################
	size(qa::QuArray, i...) = size(qa.coeffs, i...)

	getindex(qa::QuArray, i::AbstractArray) = getindex(qa.coeffs, i)
	getindex(qa::QuArray, i::Real) = getindex(qa.coeffs, i)
	getindex(qa::QuArray, i) = getindex(qa.coeffs, i)

	getindex(qa::QuArray, i...) = getindex(qa.coeffs, i...)

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

	function statevec(s::Int, fb::FiniteBasis)
		fv = quarr(zeros(Complex128, length(fb)), fb)
		fv.coeffs[s] = one(eltype(fv.coeffs))
		return fv
	end

	statevec(s::Int, lens::Int...) = statevec(s, FiniteBasis(lens))

export AbstractQuArray,
	AbstractQuVector,
	AbstractQuMatrix,
	QuArray,
	fockarr,
	statevec,
	DiracArray
