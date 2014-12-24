module QuBase
	
	#############
	# Constants #
	#############
		const lang = "\u27E8"
		const rang = "\u27E9"
		const otimes = "\u2297"
		const vdots ="\u205E"

	##################
	# Abstract Types #
	##################
		abstract AbstractStructure
		abstract AbstractQuantum{S<:AbstractStructure}
		abstract AbstractState{S<:AbstractStructure} <: AbstractQuantum{S}
		abstract AbstractOperator{S<:AbstractStructure} <: AbstractQuantum{S}
		abstract AbstractBasis{S<:AbstractStructure} <: AbstractQuantum{S}
		abstract AbstractQuArray{S<:AbstractStructure, T, N} <: AbstractArray{T,N}

	#############
	# Functions #
	#############
		structure{S}(::AbstractQuantum{S}) = S
		structure{S}(::Type{AbstractQuantum{S}}) = S
		structure(::Type{AbstractQuantum}) = AbstractStructure

		structure{S}(::Type{AbstractState{S}}) = S
		structure(::Type{AbstractState}) = AbstractStructure
		
		structure{S}(::Type{AbstractOperator{S}}) = S
		structure(::Type{AbstractOperator}) = AbstractStructure

		structure{S}(::AbstractQuArray{S}) = S
		structure{S,T,N}(::Type{AbstractQuArray{S,T,N}}) = S
		structure(::Type{AbstractQuArray}) = AbstractStructure

	######################
	# Include Statements #
	######################
		include("statelabel.jl")
		include("diracstates.jl")
		include("diracoperators.jl")

	export AbstractStructure, 
		AbstractQuantum,
		AbstractState,
		AbstractOperator,
		structure
		
end