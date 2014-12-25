import Base: 
	getindex,
	length,
	size,
	ndims,
	start,
	done,
	next,
	endof,
	last,
	first,
	collect,
	ctranspose,
	repr,
	show

#############
# FockBasis #
#############
	# A FockBasis is a wrapper around a FiniteBasis that 
	# uses precomputed values to efficiently generate 
	# a state for a given index in the basis, or
	# vice versa (an index for a given state in the basis).
	# For example:
	#
	#	julia> f=FockBasis(Ket,2,2,2)
	#	 FockBasis{Ket,AbstractStructure}(2,2,2)	
	#	
	#	julia> collect(f)
	#	 8-element Array{DiracState{Ket,AbstractStructure},1}:
	#	 | 0,0,0 ⟩
	#	 | 1,0,0 ⟩
	#	 | 0,1,0 ⟩
	#	 | 1,1,0 ⟩
	#	 | 0,0,1 ⟩
	#	 | 1,0,1 ⟩
	#	 | 0,1,1 ⟩
	#	 | 1,1,1 ⟩
	#
	# 	julia> f[6]
	# 	 | 1,0,1 ⟩
	#
	# 	julia> f[ket(1,0,1)]
	# 	 6
	# 
	# Because the states are generated rather than actually 
	# stored, one can represent very large bases without 
	# any storage overhead:
	#
	#	julia> f = FockBasis(Ket, FiniteBasis(221,135,31,42,321,3))
	#	 FockBasis{Ket,AbstractStructure}(221,135,31,42,321,3)
	#
	# 	julia> length(f)
	# 	 37407898710
	#
	# 	julia> last(f)
	# 	 | 220,134,30,41,320,2 ⟩
	#
	# 	julia> f[34234134]
	# 	 | 128,60,0,37,0,0 ⟩
	#
	# 	julia> f[ket(128,60,0,37,0,0)]
	# 	 34234134
	
	immutable FockBasis{D,S} <: AbstractDiracBasis{D,S}
		basis::FiniteBasis{S}
		denoms::(Float64...)
		FockBasis(basis, denoms, ::Type{BypassFlag}) = new(basis, denoms)
		# reverse to match cartesianmap tensor order
		# a reversal happens here *and* in precompute_denoms
		FockBasis(basis) = FockBasis{D,S}(basis, precompute_denoms(reverse(size(basis))), BypassFlag) 
	end

	FockBasis{D,S}(::Type{D}, basis::FiniteBasis{S}) = FockBasis{D,S}(basis)
	FockBasis{D}(::Type{D}, lens::Int...) = FockBasis(D, FiniteBasis(lens))

	convert{D,S}(::Type{FockBasis{D,S}}, f::FockBasis) = FockBasis{D,S}(f.basis, f.denoms, BypassFlag)

	####################
	# Helper Functions #
	####################
	# This function precomputes the 
	# denominators for each factor of 
	# the cartesian product.
	#
	# This site offers a thorough 
	# explanation of this method:
	# http://phrogz.net/lazy-cartesian-product
	#
	function precompute_denoms(lens)
		# storing as Floats avoids number precision issues for 
		# outrageously large bases
		total_divisor = prod(map(float,lens)) 
		function get_denom(i)
			total_divisor = div(total_divisor, i)
			return max(1.0, total_divisor)
		end
		return reverse(map(get_denom, lens))
	end

	ind_value(n, denom, modulus) = int(div(n, denom) % modulus)

	tuple_at_ind(f::FockBasis, i) = ntuple(ndims(f), x->ind_value(i-1, f.denoms[x], size(f.basis,x)))
	label_at_ind(f::FockBasis, i) = StateLabel(tuple_at_ind(f, i))

	######################
	# Property Functions #
	######################
	structure{D,S}(::Type{FockBasis{D,S}}) = S
	structure(::Type{FockBasis}) = AbstractStructure

	dualtype{D,S}(::Type{FockBasis{D,S}}) = D
	dualtype(::Type{FockBasis}) = DualType

	length(f::FockBasis) = length(f.basis)
	size(f::FockBasis) = size(f.basis)
	ndims(f::FockBasis) = ndims(f.basis)

	checkcoeffs(coeffs, dim, f::FockBasis) = checkcoeffs(coeffs, dim, f.basis)

	######################
	# Accessor Functions #
	######################
	getpos(f::FockBasis, s) = int(sum(map(*, gettuple(label(s)), f.denoms)))+1

	getindex{D,S}(f::FockBasis{D,S}, i) = DiracState{D,S}(label_at_ind(f, i))
	getindex(f::FockBasis, s::AbstractDiracState) = getpos(f, s)

	getindex(f::FockBasis, arr::AbstractArray) = [f[i] for i in arr]

	######################
	# Iterator Functions #
	######################
	start(::FockBasis) = 1
	done(f::FockBasis, state) = length(f) == state-1
	next(f::FockBasis, state) = f[state], state+1
	endof(f::FockBasis) = length(f)
	last(f::FockBasis) = f[length(f)]
	first(f::FockBasis) = f[1]
	collect(f::FockBasis) = f[1:end]
	collectlabels(f::FockBasis) = [label_at_ind(f, i) for i=1:length(f)]

	##########################
	# Mathematical Functions #
	##########################
	tensor{D}(a::FockBasis{D}, b::FockBasis{D}) = FockBasis(D, tensor(a.basis, b.basis))
	ctranspose{D,S}(f::FockBasis{D,S}) = convert(FockBasis{D',S}, f)

	######################
	# Printing Functions #
	######################
	repr(f::FockBasis) = "$(typeof(f))$(size(f.basis))"
	show(io::IO, f::FockBasis) = print(io, repr(f))

export FockBasis,
	structure,
	dualtype,
	getpos