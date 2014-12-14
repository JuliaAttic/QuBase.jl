module QuBase

	abstract AbstractSpace
	abstract AbstractState{S<:AbstractSpace}
	abstract AbstractOperator{S<:AbstractSpace}
	abstract AbstractQuArray{S<:AbstractSpace, T, N} <: AbstractArray{T,N}

	export AbstractSpace, 
		AbstractState,
		AbstractOperator,
		AbstractQuArray
		
end