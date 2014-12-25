import Base:
	convert,
	copy,
	size,
	length,
	in,
	getindex,
	filter,
	map,
	ctranspose,
	vcat,
	setdiff,
	kron,
	summary,
	show

####################
# Helper Functions #
####################
	combine_hashes(hash1, hash2) = hash1 $ hash2

##############
# DiracBasis #
##############
	immutable DiracBasis{D<:DualType, S<:AbstractStructure} <: AbstractDiracBasis{D,S}
		labels::Vector{StateLabel}
		labelmap::Dict{StateLabel, Int}
		labels_hash::Uint64
		
		function DiracBasis(labels::AbstractArray{StateLabel}, labelmap::Dict{StateLabel, Int}, labels_hash::Uint64)
			return new(labels, labelmap, labels_hash)
		end 

		function DiracBasis(labels::AbstractArray{StateLabel}, ::Type{BypassFlag})
			labelmap = Dict{StateLabel,Int}()
			sizehint(labelmap, length(labels))
			
			labelmap[labels[1]] = 1
			labels_hash = hash(labels[1])
			for i=2:length(labels)
				labels_hash = combine_hashes(labels_hash, hash(labels[i]))
				labelmap[labels[i]] = i
			end

			return DiracBasis{D,S}(labels, labelmap, labels_hash)
		end
		
		DiracBasis(labels::AbstractArray{StateLabel}) = DiracBasis{D,S}(unique(labels), BypassFlag)
	end

	DiracBasis{D,S}(::Type{D}, ::Type{S}, labels::AbstractArray{StateLabel}, flag::Type{BypassFlag}) = DiracBasis{D,S}(labels, flag)
	DiracBasis{D,S}(::Type{D}, ::Type{S}, labels::AbstractArray{StateLabel}) = DiracBasis{D,S}(labels)
	DiracBasis{D,S}(::Type{D}, ::Type{S}, labels::(AbstractArray{StateLabel}...,)) = DiracBasis(D, S, cart_prod(labels))
	DiracBasis{D,S}(::Type{D}, ::Type{S}, arr::AbstractArray) = DiracBasis(D, S, map(StateLabel, arr))
	DiracBasis{D,S}(::Type{D}, ::Type{S}, arr::AbstractArray...) = DiracBasis(D, S, map(arr->map(StateLabel, arr), arr))
	DiracBasis{D,S}(::Type{D}, ::Type{S}, labels...) = DiracBasis(D, S, collect(labels))
	 
	DiracBasis{D,S}(arr::AbstractArray{AbstractDiracState{D,S}}) = DiracBasis(D, S, arr)
	DiracBasis{D}(arr::AbstractArray{AbstractDiracState{D}}) = DiracBasis(D, reduce(typejoin, map(structure, arr)), arr)
	DiracBasis(arr::AbstractArray{AbstractDiracState}) = DiracBasis(Ket, reduce(typejoin, map(structure, arr)), arr)
	DiracBasis(states::AbstractDiracState...) = DiracBasis(collect(states))

	DiracBasis(labels...) = DiracBasis(Ket, AbstractStructure, labels...)

	convert{D,S}(::Type{DiracBasis{D,S}}, b::DiracBasis) = DiracBasis{D,S}(b.labels, b.labelmap, b.labels_hash)
	copy{D,S}(b::DiracBasis{D,S}) = DiracBasis{D,S}(copy(b.labels), copy(b.labelmap), copy(b.labels_hash))

	######################
	# Property Functions #
	######################
	structure{D,S}(::Type{DiracBasis{D,S}}) = S
	structure(::Type{DiracBasis}) = AbstractStructure

	dualtype{D,S}(::Type{DiracBasis{D,S}}) = D
	dualtype(::Type{DiracBasis}) = DualType

	labels(b::DiracBasis) = b.labels
	samelabels(a::DiracBasis, b::DiracBasis) = a.labels_hash == b.labels_hash
	nfactors(B::DiracBasis) = length(first(B.labels))
	size(B::DiracBasis) = size(B.labels)
	length(B::DiracBasis) = length(B.labels)

	checkcoeffs(coeffs, dim, basis::DiracBasis) = size(coeffs, dim) == length(basis) 

	########################
	# Array-like Functions #
	########################
	getpos(B::DiracBasis, label::StateLabel) = B.labelmap[label]
	in(label::StateLabel, b::AbstractDiracBasis) = haskey(b.labelmap, label)

	getindex(B::DiracBasis, s::AbstractDiracState) = getpos(B, label(s))
	generic_getindex{D,S}(B::DiracBasis{D,S}, i) = DiracState{D, S}(getindex(labels(B), i))

	generic_collect{D,S}(B::DiracBasis{D,S}, arr) = DiracBasis{D,S}(B.labels[arr], BypassFlag)

	getindex(B::DiracBasis, arr::AbstractArray) = generic_collect(B, arr)
	getindex(B::DiracBasis, i) = generic_getindex(B, i)

	filter{D,S}(f::Function, B::DiracBasis{D,S}) = DiracBasis(D, S, filter(f, B.labels), BypassFlag)
	map{D,S}(f::Function, B::DiracBasis{D,S}) = DiracBasis(D, S, map(f, B.labels))

	filter{D,S}(f::Function, B::FockBasis{D,S}) = DiracBasis(D, S, filter(f, collectlabels(B)), BypassFlag)
	map{D,S}(f::Function, B::FockBasis{D,S}) = DiracBasis(D, S, map(f, collectlabels(B)))
	
	xsubspace(B::AbstractDiracBasis, x::Int) = filter(s->sum(gettuple(s))==x, B)

	#####################
	# Joining Functions #
	#####################
	function append_label!(map::Dict, s::StateLabel)
		map[s] = length(map)+1
		return map
	end

	function prepend_label!(map::Dict, s::StateLabel)
		map[s] = 1
		for (k, v) in map
			map[k] = v+1
		end
		return map
	end

	non_unique_vcat{D,S}(s::StateLabel, b::DiracBasis{D,S}) = DiracBasis{D,S}(vcat(b.labels, s), prepend_label!(copy(b.labelmap), s, ), combine_hashes(hash(s), b.labels_hash))
	non_unique_vcat{D,S}(b::DiracBasis{D,S}, s::StateLabel) = DiracBasis{D,S}(vcat(b.labels, s), append_label!(copy(b.labelmap), s), combine_hashes(b.labels_hash, hash(s)))

	vcat{D,S}(s::StateLabel, b::DiracBasis{D,S}) = in(s, b) ? non_unique_vcat(s, b) : b
	vcat{D,S}(b::DiracBasis{D,S}, s::StateLabel) = in(s, b) ? non_unique_vcat(s, b) : b

	function vcat{D,A,B}(a::DiracBasis{D,A}, b::DiracBasis{D,B}) 
		res = convert(DiracBasis{D, typejoin(A,B)}, a)
		for label in b.labels
			res = vcat(res, label)
		end
		return res
	end

	vcat{D,A,B}(s::AbstractDiracState{D,A}, b::AbstractDiracBasis{D,B}) = vcat(label(s), convert(DiracBasis{D, typejoin(A,B)}, b))
	vcat{D,A,B}(b::AbstractDiracBasis{D,A}, s::AbstractDiracState{D,B}) = vcat(convert(DiracBasis{D, typejoin(A,B)}, b), label(s))

	setdiff{D, A<:AbstractStructure,B<:AbstractStructure}(a::AbstractDiracBasis{D,A}, b::AbstractDiracBasis{D,B}) = DiracBasis{D,typejoin(A,B)}(setdiff(a.labels, b.labels), BypassFlag)

	##########################
	# Mathematical Functions #
	##########################
	function cart_prod{A<:AbstractVector{StateLabel}}(labels::(A...,))
		lens = map(length, labels)
		arr = Array(StateLabel, prod(lens))
		index = 1

		function set_ind!(inds...)
			arr[index] = combine(map(getindex, labels, inds))
			index += 1
		end
		
		cartesianmap(set_ind!, lens)
		return arr
	end

	kron{D}(bases::AbstractDiracBasis{D}...) = DiracBasis{D,reduce(typejoin, map(structure, bases))}(cart_prod(map(labels, bases)), BypassFlag)	
	kron{D}(basis::AbstractDiracBasis{D}, s::AbstractDiracState{D}) = kron(basis, DiracBasis(s))	
	kron{D}(s::AbstractDiracState{D}, basis::AbstractDiracBasis{D}) = kron(DiracBasis(s), basis)	

	dir_prod(a::Vector{StateLabel}, b::Vector{StateLabel}) =  length(a)==length(b) ? StateLabel[flatten(a[i], b[i]) for i=1:length(a)] : error("Arrays must be same length")
	directprod{D}(bases::AbstractDiracBasis{D}...) = DiracBasis{D, reduce(typejoin, map(structure, bases))}(reduce(dir_prod, map(labels, bases)), BypassFlag)

	factorize{D,S}(basis::DiracBasis{D,S}) = [DiracBasis{D,S}([StateLabel(basis[i][j]) for i=1:length(basis)]) for j=1:nfactors(basis)]

	ctranspose{S<:AbstractStructure}(B::DiracBasis{Ket,S}) = convert(DiracBasis{Bra,S}, B)
	ctranspose{S<:AbstractStructure}(B::DiracBasis{Bra,S}) = convert(DiracBasis{Ket,S}, B)

	######################
	# Printing Functions #
	######################
	summary(basis::DiracBasis) = "$(typeof(basis)) with $(length(basis)) states:"
	function show(io::IO, basis::DiracBasis)
		print(io, summary(basis))
		
		if length(basis) > 25
			
			println(io)
			for i=1:10
				println(io, basis[i])
			end

			print(io, vdots) # vdots is repo-wide const
			
			start_ind = length(basis)-10
		else
			start_ind = 1
		end

		for i=start_ind:length(basis)
			println(io)
			print(io, basis[i])
		end
	end

export DiracBasis,
	structure,
	dualtype,
	nfactors,
	factorize,
	getpos,
	samelabels,
	directprod,
	xsubspace