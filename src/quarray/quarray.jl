import Base: 
	size,
	getindex

abstract AbstractQuArray{B<:(AbstractBasis...), T, N} #<: AbstractArray{T,N}

typealias AbstractQuVector{B<:(AbstractBasis), T} AbstractQuArray{B, T}
typealias AbstractQuMatrix{B<:(AbstractBasis, AbstractBasis), T} AbstractQuArray{B, T}

###########
# QuArray # 				
###########

	checkcoeffs(coeffs, bases::(AbstractBasis...)) = reduce(&&, [checkcoeffs(coeffs, i, bases[i]) for i=1:ndims(coeffs)])

	type QuArray{B<:(AbstractBasis...), T, N, A} <: AbstractQuArray{B, T, N}
	    coeffs::A
	    bases::B
	    function QuArray(coeffs::AbstractArray{T, N}, bases::B) 
	    	if N == length(bases) && checkcoeffs(coeffs, bases) 
	    		new(coeffs, bases)
	    	else 
	    		error("Coefficient array does not conform to input bases")
	    	end
	    end
	end

	QuArray{T,N}(coeffs::AbstractArray{T,N}, bases::Tuple) = QuArray{typeof(bases), T, N, typeof(coeffs)}(coeffs, bases)
	QuArray(coeffs, bases...) = QuArray(coeffs, bases)

	size(qa::QuArray, i...) = size(qa.coeffs, i...)
	length(qa::QuArray) = length(qa.coeffs)
	getindex(qa::QuArray, i...) = getindex(qa.coeffs, i...)


# # standard operations with QuArrays
# *{B<:AbstractBasis}(qa1::QuArray{B}, qa2::QuArray{B}) = QuArray(qa1.coeffs*qa2.coeffs, qa2.bases)

# # like with Arrays we can alias QuArray for different purposes

# # examples for convenience constructors
# function statevec(fb::FiniteBasis, s::Int)
# 	fsv = QuArray(zeros(Complex128, fb.nb), fb)
# 	fsv.coeffs[s] = one(eltype(fsv.coeffs))
# 	return fsv
# end

# statevec(nb::Int, s::Int) = statevec(FiniteBasis(nb), s)

# function creationop(fb::FiniteBasis)
# 	nb = fb.nb
# 	co = QuArray(sparse( [2:nb], [1:nb-1], complex(sqrt(linspace( 1, nb-1, nb-1))), nb, nb ),
# 	fb)
# 	return co
# end

# creationop(nb::Int) = creationop(FiniteBasis(nb))

# include("diracarray.jl")

export AbstractQuArray,
	AbstractQuVector,
	AbstractQuMatrix,
	QuArray
