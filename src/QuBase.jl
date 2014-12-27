module QuBase
	
	import Base: kron
		
	####################
	# String Constants #
	####################
		const lang = "\u27E8"
		const rang = "\u27E9"
		const otimes = "\u2297"
		const vdots ="\u205E"

	##################
	# Abstract Types #
	##################
		abstract AbstractStructure
		abstract AbstractQuantum{S<:AbstractStructure}

		# Various constructor methods in this repo allow an argument 
		# of type Type{BypassFlag} to be passed in in order to 
		# circumvent value precalculation/checking. This is useful for
		# conversion methods and the like, where you know the input 
		# has already been vetted elsewhere. Don't use this unless
		# you're sure of what you're doing, and don't export this.
		abstract BypassFlag

	######################
	# Include Statements #
	######################
		include("dirac/dirac.jl")
		include("bases/bases.jl")
		include("arrays/quarray.jl")
	
	#############
	# Functions #
	#############
		structure{S}(::AbstractQuantum{S}) = S
		structure{S}(::Type{AbstractQuantum{S}}) = S
		structure(::Type{AbstractQuantum}) = AbstractStructure

		# an n-arity form of the tensor
		# product, reduction is done via
		# the binary definition of tensor()
		# defined in the files included above. 
		tensor(s...) = reduce(tensor, s)

		# For the sake of convenience, kron() 
		# is defined to be equivalent to 
		# tensor() for quantum objects
		kron(a::AbstractQuantum, b::AbstractQuantum) = tensor(a, b)

	export AbstractStructure, 
		AbstractQuantum,
		structure,
		tensor	
end
	