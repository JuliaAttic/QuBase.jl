abstract AbstractBasis{S<:AbstractStructure} <: AbstractQuantum{S}
abstract AbstractLabelBasis{S<:AbstractStructure} <: AbstractBasis{S}

structure{S}(::Type{AbstractBasis{S}}) = S
structure(::Type{AbstractBasis}) = AbstractStructure

structure{S}(::Type{AbstractLabelBasis{S}}) = S
structure(::Type{AbstractLabelBasis}) = AbstractStructure

# All basis types should implement: 
#
# 	checkcoeffs(coeffs, dim, basis::NewBasisType) -> Bool
#
# ...which is then used by QuArray to ensure that 
# the coefficient array is not malformed with respect
# to the input bases. The second argument, `dim`, specifies the 
# the dimension of the coefficient array which 
# corresponds to the provided basis.

include("finitebasis.jl")
include("fockbasis.jl")
