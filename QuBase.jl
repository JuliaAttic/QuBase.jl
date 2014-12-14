module QuBase

	abstract AbstractSpace
	abstract AbstractQuantum{S<:AbstractSpace}
	abstract AbstractState{S<:AbstractSpace} <: AbstractQuantum{S}
	abstract AbstractOperator{S<:AbstractSpace} <: AbstractQuantum{S}

	space{S}(::AbstractQuantum{S}) = S
	space{S}(::Type{AbstractQuantum{S}}) = S
	space(::Type{AbstractQuantum}) = AbstractSpace

	export AbstractSpace, 
		AbstractQuantum,
		AbstractState,
		AbstractOperator,
		space
		
end