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
		abstract QuantumScalar <: Number

	#############
	# Functions #
	#############
		structure{S}(::AbstractQuantum{S}) = S

		for T=(:AbstractQuantum, :AbstractState, :AbstractOperator, :AbstractBasis)
			@eval begin
				structure{S}(::Type{($T){S}}) = S
				structure(::Type{($T)}) = AbstractStructure
			end
		end

	######################
	# Include Statements #
	######################
		include("statelabel.jl")
		include("diracstates.jl")
		include("diracoperators.jl")
		include("scalar.jl")
		include("quarray.jl")

	export AbstractStructure, 
		AbstractQuantum,
		AbstractState,
		AbstractOperator,
		AbstractBasis,
		QuantumScalar,
		structure
		
end
	