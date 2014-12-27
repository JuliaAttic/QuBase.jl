# Dirac arrays are wrappers around QuArrays which generate
# ScaledStates and ScaledOperators as elements. These 
# quantum elements are formed by multiplying together 
# the underlying QuArray's basis states with the associated
# coefficients. 
#
# For example, the nth element of DiracVector{Ket, S} is 
# the nth coefficient times a DiracState{Ket, S} whose label
# is the nth label in the basis. 
#
# Likewise, the (ith, jth) element of DiracMatrix{Ket, S} is 
# the (ith, jth) coefficient times a DiracOperator{S} whose 
# ket label is the ith label in the row basis, and whose 
# bra label is the jth label in the column basis.
#
# This kind of structure means we can do two cool things:
# 
# 1) If operations are being done within a single basis, we can 
# bypass the performance cost of using states/basis and 
# just operate on the coefficient arrays. If we're doing 
# mixed basis operations, however, we can perform operations
# utilizing the bras and kets themselves. 
#
# 2) Storing a basis of labels allows for label-based methods 
# of analysis. An easy example is arbitrary selection/extraction 
# of subspaces using methods like `filter`. (not yet implemented
# here, but examples of which can be found in the old QuDirac repo).

abstract DiracArray{B, T<:AbstractDirac, N} <: AbstractQuArray{B, T, N}

###############
# DiracVector #
###############
	type DiracVector{D, S<:AbstractStructure, T, B<:AbstractLabelBasis} <: DiracArray{(B,), ScaledState{D, S, T}, 1}
		arr::QuArray{(B,), T}
		DiracVector(arr::QuArray{(AbstractLabelBasis{S},), T}) = new(arr)
	end

###############
# DiracMatrix #
###############
	type DiracMatrix{S<:AbstractStructure, 
					 T, 
					 R<:AbstractLabelBasis, 
					 C<:AbstractLabelBasis} <: DiracArray{(R,C), ScaledOperator{S, T}, 2}
		arr::QuArray{(R,C), T}
		DiracMatrix(arr::QuArray{(AbstractLabelBasis{S},AbstractLabelBasis{S}), T, N, A}) = new(arr)
	end

export DiracArray,
	DiracVector,
	DiracMatrix