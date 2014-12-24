import Base: *

abstract AbstractQuArray{B<:(AbstractBasis...), T, N} <: AbstractArray{T,N}

# basic basis consisting of a finite number of states
# should the number of basis states be a type parameter?
immutable FiniteBasis <: AbstractBasis
	nb::Int		# number of vectors in basis
end

abstract AbstractProductBasis <: AbstractBasis

# we also need product states, for example if one has different subsystems
# one possiblity would be to collect the basis for each of those
type FinitePBasis <: AbstractProductBasis
	bvec::Vector{FiniteBasis}
end

# the work horse for all computations
type QuArray{B<:(AbstractBasis...), T, N, A<:AbstractArray} <: AbstractQuArray{B, T, N}
    coeffs::A
    bases::B
    function QuArray(coeffs::AbstractArray{T, N}, basis::(AbstractBasis...))
        # size/dimension checking should go here
        new(coeffs, basis)
    end
end

QuArray{T,N}(coeffs::AbstractArray{T,N}, bases::AbstractBasis...) = QuArray{typeof(bases), T, N, typeof(coeffs)}(coeffs, bases)
QuArray(c::AbstractArray) = QuArray(c, FiniteBasis(size(c,1)))

# standard operations with QuArrays
*{B<:AbstractBasis}(qa1::QuArray{B}, qa2::QuArray{B}) = QuArray(qa1.elem*qa2.elem, qa2.basis)

# like with Arrays we can alias QuArray for different purposes
typealias QuStateVector{S, B<:AbstractBasis, T} QuArray{B, T, 1}
typealias QuOperator{S, B<:AbstractBasis, T} QuArray{B, T, 2}

# examples for convenience constructors
function statevec(fb::FiniteBasis, s::Int)
	fsv = QuArray(zeros(Complex128, fb.nb), fb)
	fsv.elem[s] = one(eltype(fsv.elem))
	return fsv
end

statevec(nb::Int, s::Int) = statevec(FiniteBasis(nb), s)

function creationop(fb::FiniteBasis)
	nb = fb.nb
	co = QuArray(sparse( [2:nb], [1:nb-1], complex(sqrt(linspace( 1, nb-1, nb-1))), nb, nb ),
	fb)
	return co
end

creationop(nb::Int) = creationop(FiniteBasis(nb))

export AbstractQuArray,
	FiniteBasis,
	FinitePBasis,
	QuArray,
	QuStateVector,
	QuOperator,
	statevec,
	creationop
