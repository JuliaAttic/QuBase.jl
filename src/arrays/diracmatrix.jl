type DiracMatrix{S<:AbstractStructure, T, N, A} <: AbstractQuArray{(FockBasis{S},FockBasis{S}), ScaledOperator{S, T}, N}
	arr::QuArray{(FockBasis{S},FockBasis{S}), T, N, A}
end