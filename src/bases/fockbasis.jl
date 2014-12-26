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
	show,
	in

#############
# FockBasis #
#############
	# A FockBasis is a wrapper around a FiniteBasis that 
	# uses precomputed values to efficiently generate 
	# StateLabels for given indices in the basis, or
	# vice versa (an index for a given state in the basis).
	# For example:
	#
	#	julia> f=FockBasis(2,2,2)
	#	 FockBasis{AbstractStructure}(2,2,2)	
	#	
	#	julia> collect(f)
	#	 8-element Array{StateLabel,1}:
	#	 StateLabel(0,0,0)
	#	 StateLabel(1,0,0)
	#	 StateLabel(0,1,0)
	#	 StateLabel(1,1,0)
	#	 StateLabel(0,0,1)
	#	 StateLabel(1,0,1)
	#	 StateLabel(0,1,1)
	#	 StateLabel(1,1,1)
	#
	# 	julia> f[6]
	# 	 StateLabel(1,0,1)
	#
	# 	julia> f[StateLabel(1,0,1)]
	# 	 6
	# 
	# Because the labels are generated rather than actually 
	# stored, one can represent very large bases without 
	# any storage overhead:
	#
	#	julia> f = FockBasis(221,135,31,42,321,3)
	#	 FockBasis{AbstractStructure}(221,135,31,42,321,3)
	#
	# 	julia> length(f)
	# 	 37407898710
	#
	# 	julia> last(f)
	# 	 StateLabel(220,134,30,41,320,2)
	#
	# 	julia> f[34234134]
	# 	 StateLabel(128,60,0,37,0,0)
	#
	# 	julia> f[StateLabel(128,60,0,37,0,0)]
	# 	 34234134
	
	immutable FockBasis{S<:AbstractStructure} <: AbstractBasis{S}
		basis::FiniteBasis{S}
		denoms::(Float64...)
		FockBasis(basis, denoms, ::Type{BypassFlag}) = new(basis, denoms)
		# reverse is done to match cartesianmap order
		FockBasis(basis::FiniteBasis{S}) = FockBasis{S}(basis, precompute_denoms(reverse(size(basis))), BypassFlag) 
		FockBasis(lens::(Int...)) = FockBasis{S}(FiniteBasis{S}(lens))
		FockBasis(lens::Int...) = FockBasis{S}(FiniteBasis{S}(lens))
	end

	FockBasis{S}(basis::FiniteBasis{S}) = FockBasis{S}(basis)
	FockBasis(lens::(Int...)) = FockBasis(FiniteBasis(lens)) 
	FockBasis(lens::Int...) = FockBasis(FiniteBasis(lens))

	convert{S}(::Type{FockBasis{S}}, f::FockBasis) = FockBasis{S}(convert(FiniteBasis{S}, f.basis), f.denoms, BypassFlag)

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
		# reverse is done to match cartesianmap order
		return reverse(map(get_denom, lens))
	end

	ind_value(n, denom, modulus) = int(div(n, denom) % modulus)

	tuple_at_ind(f::FockBasis, i) = ntuple(ndims(f), x->ind_value(i-1, f.denoms[x], size(f.basis,x)))

	######################
	# Property Functions #
	######################
	structure{S}(::Type{FockBasis{S}}) = S
	structure(::Type{FockBasis}) = AbstractStructure

	length(f::FockBasis) = length(f.basis)
	size(f::FockBasis) = size(f.basis)
	ndims(f::FockBasis) = ndims(f.basis)

	checkcoeffs(coeffs, dim, f::FockBasis) = checkcoeffs(coeffs, dim, f.basis)

	######################
	# Accessor Functions #
	######################
	checkrange(x,y) = 0 <= x < y
	in(s, f::FockBasis) = reduce(&, map(checkrange, s, size(f.basis)))
	getpos(f::FockBasis, s) = int(sum(map(*, s, f.denoms)))+1 

	getindex(f::FockBasis, i) = StateLabel(tuple_at_ind(f, i))
	getindex(f::FockBasis, s::StateLabel) = s in f ? getpos(f, s) : error("StateLabel not found: $s")
	getindex(f::FockBasis, arr::AbstractArray) = [f[i] for i in arr]

	getstate{D<:DualType, S}(f::FockBasis{S}, i, ::Type{D}=Ket) = DiracState{D, S}(f[i])
	getstate{D<:DualType}(f::FockBasis, arr::AbstractArray, ::Type{D}=Ket) = [getstate(f, i, D) for i in arr]

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

	##########################
	# Mathematical Functions #
	##########################
	tensor(a::FockBasis, b::FockBasis) = FockBasis(tensor(a.basis, b.basis))

	######################
	# Printing Functions #
	######################
	repr(f::FockBasis) = "$(typeof(f))$(size(f.basis))"
	show(io::IO, f::FockBasis) = print(io, repr(f))

export FockBasis,
	structure,
	getstate,
	tensor