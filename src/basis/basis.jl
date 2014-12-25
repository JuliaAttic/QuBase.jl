abstract AbstractBasis{S<:AbstractStructure} <: AbstractQuantum{S}
abstract AbstractDiracBasis{D<:DualType, S<:AbstractStructure} <: AbstractBasis{S}

structure{S}(::Type{AbstractBasis{S}}) = S
structure(::Type{AbstractBasis}) = AbstractStructure

structure{S}(::Type{AbstractDiracBasis{S}}) = S
structure(::Type{AbstractDiracBasis}) = AbstractStructure

dualtype{D}(::AbstractDiracBasis{D}) = D
dualtype{D,S}(::Type{AbstractDiracBasis{D,S}}) = D
dualtype(::Type{AbstractDiracBasis}) = DualType

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
include("diracbasis.jl")
