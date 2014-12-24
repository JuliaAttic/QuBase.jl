module QuBase

	import Base: *

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

	###########################################################################
	abstract AbstractBasis

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


	abstract AbstractQuArray{B<:AbstractBasis,N}

	# the work horse for all computations; note that there is only
	# one basis specfied, maybe this should/could be generalized?
	type QuArray{B<:AbstractBasis,N,CT,T} <: AbstractQuArray{B,N}
		elem::CT	# container for coefficients
		basis::B	# basis associated with coefficients

		function QuArray(c::AbstractArray{T,N}, b::B)
			# should check that the size of c and b are compatible
			new(c,b)
		end
	end

	# outer constructor
	QuArray{B<:AbstractBasis,N,T}(c::AbstractArray{T,N}, b::B) = QuArray{B,N,typeof(c),T}(c, b)
	QuArray(c::AbstractArray) = QuArray(c, FiniteBasis(size(c,1)))

	# standard operations with QuArrays
	*{B<:AbstractBasis}(qa1::QuArray{B}, qa2::QuArray{B}) = QuArray(qa1.elem*qa2.elem, qa2.basis)

	# like with Arrays we can alias QuArray for different purposes
	typealias QuStateVector{B<:AbstractBasis} QuArray{B,1}
	typealias QuOperator{B<:AbstractBasis} QuArray{B,2}


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
end
