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
		abstract QuantumScalar <: Number

		# Various constructor methods in this repo allow an argument 
		# of type Type{BypassFlag} to be passed in in order to 
		# circumvent value precalculation/checking; this is useful for
		# conversion methods and the like. Don't use this unless
		# you're sure of what you're doing, and DEFINITELY 
		# don't export this.
		abstract BypassFlag

	#############
	# Functions #
	#############
		structure{S}(::AbstractQuantum{S}) = S

		for T=(:AbstractQuantum, :AbstractState, :AbstractOperator)
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
		include("basis/basis.jl")
		# include("quarray/quarray.jl")

	export AbstractStructure, 
		AbstractQuantum,
		AbstractState,
		AbstractOperator,
		AbstractBasis,
		QuantumScalar,
		structure
		
end
	