type DiracVector{D, S<:AbstractStructure, T, N, A} <: AbstractQuArray{(FockBasis{S},), ScaledState{D, S, T}, N}
	arr::QuArray{(FockBasis{S},), T, N, A}
end
