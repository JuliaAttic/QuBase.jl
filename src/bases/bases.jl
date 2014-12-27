####################
# Type Definitions #
####################
	abstract AbstractBasis{S<:AbstractStructure} <: AbstractQuantum{S}
	abstract AbstractLabelBasis{S<:AbstractStructure} <: AbstractBasis{S}

######################
# Include Statements #
######################
	include("finitebasis.jl")
	include("fockbasis.jl")
	include("labelbasis.jl")

#############
# Functions #
#############
	# Note: All B<:AbstractBasis types should implement the following: 
	#
	# 	checkcoeffs(coeffs, dim, basis::B) -> Bool
	#
	# ...used by QuArray to ensure that the coefficient array is 
	# not malformed with respect to the input bases. The second 
	# argument, `dim`, specifies the dimension of the coefficient 
	# array which corresponds to the provided basis.
	#
	# All B<:AbstractLabelBasis types should implement the following:
	# 
	#	getindex(basis::B, i) -> the StateLabel at index `i`
	# 	getindex(basis::B, label::StateLabel) -> the index at which `label` resides in `basis`
	#	in(label::StateLabel, basis::B) -> checks if the given `label` is an element of `basis`
	#   nfactors(basis::B) -> the number of factor bases for `basis`; this is the same as the 
	#						  length of a StateLabel in the basis, or the number of particles
	#						  present in this basis.
	#	labelvec(basis::B) -> returns a Vector{StateLabel} of the labels in this basis.
	#	samelabels(a::B, b::B) -> returns true if the labels of `a` are the same as 
	#							  those in `b`, and in the same order. This
	# 							  should be implmented to run in constant time/low-input 
	#                             linear time if at all possible, as this function is used 
	#							  by DiracArrays to check whether or not an array 
	# 							  operation can be performed solely with coefficients or 
	#							  requires the use of bases.  

	checkcoeffs(coeffs::AbstractArray, dim::Int, basis::AbstractBasis) = error("checkcoeffs(coeffs, dim, ::$B) must be defined!")
	structure{S}(::Type{AbstractBasis{S}}) = S
	structure(::Type{AbstractBasis}) = AbstractStructure
	structure{S}(::Type{AbstractLabelBasis{S}}) = S
	structure(::Type{AbstractLabelBasis}) = AbstractStructure

	function getstate{D<:DualType, S}(basis::AbstractLabelBasis{S}, 
		                     		  i, 
		                     		  ::Type{D}=Ket)
		return DiracState{D,S}(basis[i])
	end

	function getstate{D<:DualType}(basis::AbstractLabelBasis, 
		                  		   arr::AbstractArray, 
		                  		   ::Type{D}=Ket)
		return [getstate(basis, i, D) for i in arr]
	end

	filter{S}(f::Function, basis::AbstractLabelBasis{S}) = LabelBasis{S}(filter(f, labelvec(basis)), BypassFlag)
	map{S}(f::Function, basis::AbstractLabelBasis{S}) = LabelBasis{S}(map(f, labelvec(basis)))
	
	xsubspace(basis::AbstractLabelBasis, x::Int) = filter(s->sum(s)==x, basis)

export checkcoeffs,
	structure,
	xsubspace


