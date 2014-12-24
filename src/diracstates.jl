import Base:
	getindex,
	length,
	show,
	repr,
	convert,
	ctranspose,
	kron,
	Ac_mul_B, #TODO: Define other members of the A_mul_B family as necessary 
	+,
	-

################################
# Abstract types and Functions #
################################
	# An AbstractDiracState is a type that
	# represents a state formulated in
	# Dirac notation - an abstract vector (bra/ket)
	# optionally multiplied by a scalar value. 
	#
	# It has two type parameters. The first is a subtype
	# DualType that specifies whether the state is 
	# an element of ket-space or bra-space. The second 
	# is a subtype of AbstractStructure, which provides
	# type information regarding the structure of the basis
	# that the state belongs to.
	#
	# Concrete subtypes of AbstractDiracState 
	# implement the following methods:
	#
	# coeff(s::AbstractDiracState) -> returns the scalar coefficient of the state
	# state(s::AbstractDiracState) -> returns only the state, without a scalar coefficient, as a DiracState
	# label(s::AbstractDiracState) -> returns the state's label as a StateLabel
	# repr(s::AbstractDiracState) -> returns a string representation of the state
	# coefftype(s::AbstractDiracState), 
	# coefftype(::Type{AbstractDiracState}) -> returns the coefficient type for a given state/type
	# structure(::Type{AbstractDiracState}) -> returns the state type's structure type parameter
	# dualtype(::Type{AbstractDiracState}) -> returns the state type's dual type parameter

	abstract DualType
	abstract Ket <: DualType
	abstract Bra <: DualType

	abstract AbstractDiracState{D<:DualType, S<:AbstractStructure} <: AbstractState{S}

	typealias AbstractDiracKet{S<:AbstractStructure} AbstractDiracState{Ket, S}
	typealias AbstractDiracBra{S<:AbstractStructure} AbstractDiracState{Bra, S}

	structure{D,S}(::Type{AbstractDiracState{D,S}}) = S
	structure(::Type{AbstractDiracState}) = AbstractStructure

	dualtype{D,S}(::AbstractDiracState{D,S}) = D
	dualtype{D,S}(::Type{AbstractDiracState{D,S}}) = D
	dualtype(::Type{AbstractDiracState}) = DualType

	getindex(s::AbstractDiracState, i) = getindex(label(s), i)
	length(s::AbstractDiracState) = length(label(s))
	show(io::IO, s::AbstractDiracState) = print(io, repr(s))

	# This comparison method is useful to define
	# for all labeled concrete subtypes of AbstractQuantum
	samelabels(a::AbstractDiracState, b::AbstractDiracState) = label(a) == label(b)

##############
# DiracState #
##############
	# A DiracState is the type representation of
	# an unscaled abstract vector, formulated in 
	# Dirac notation as a bra or ket.
	immutable DiracState{D<:DualType, S<:AbstractStructure} <: AbstractDiracState{D, S}
		label::StateLabel
	end

	typealias DiracKet{S} DiracState{Ket, S}
	typealias DiracBra{S} DiracState{Bra, S}

	convert{D,S}(::Type{DiracState{D,S}}, s::AbstractDiracState) = DiracState{D,S}(label(s))

	############################
	# Convenience Constructors #
	############################
	ket(label::StateLabel, S=AbstractStructure) = DiracKet{S}(label)
	ket(tup::Tuple, S=AbstractStructure) = ket(StateLabel(tup),S)
	ket(labels...) = ket(labels)
	bra(label::StateLabel, S=AbstractStructure) = DiracBra{S}(label)
	bra(tup::Tuple, S=AbstractStructure) = bra(StateLabel(tup),S)
	bra(labels...) = bra(labels)

	################################
	# AbstractDiracState Functions #
	################################
	# We somewhat arbitrarily 
	# default to Int as the 
	# coefficient type. Then 
	# coeff can just return the multiplicative
	# identity of that chosen type - 1::Int
	coefftype(s::DiracState) = Int
	coefftype{S<:DiracState}(::Type{S}) = Int

	dualtype(::Type{DiracState}) = DualType
	dualtype{D,S}(::Type{DiracState{D,S}}) = D

	structure(::Type{DiracState}) = AbstractStructure
	structure{D,S}(::Type{DiracState{D,S}}) = S

	coeff(s::DiracState) = one(coefftype(s))
	state(s::DiracState) = s
	label(s::DiracState) = s.label

	######################
	# Printing Functions #
	######################
	labelstr(s::DiracState) = "$(labelstr(label(s)))"
	repr(k::DiracKet) = "| $(labelstr(k)) $rang"
	repr(b::DiracBra) = "$lang $(labelstr(b)) |"

	###########################
	# Mathematical Operations #
	###########################
	ctranspose{S}(k::DiracKet{S}) = DiracBra{S}(label(k))
	ctranspose{S}(b::DiracBra{S}) = DiracKet{S}(label(b))

###############
# ScaledState #
###############
	# A ScaledState is the type representation of
	# a DiracState multiplied by a scalar. 
	immutable ScaledState{D<:DualType, S<:AbstractStructure, T} <: AbstractDiracState{D, S}
		coeff::T
		state::DiracState{D,S}
	end
		
	convert{D,S,T}(::Type{ScaledState{D,S,T}}, s::AbstractDiracState) = ScaledState(convert(T, coeff(s)), convert(DiracState{D,S}, state(s)))

	typealias ScaledKet{S, T} ScaledState{Ket, S, T}
	typealias ScaledBra{S, T} ScaledState{Bra, S, T}

	######################
	# Accessor Functions #
	######################
	coefftype{D,S,T}(::ScaledState{D,S,T}) = T
	coefftype(::Type{ScaledState}) = Any
	coefftype{D,S,T}(::Type{ScaledState{D,S,T}}) = T

	dualtype(::Type{ScaledState}) = DualType
	dualtype{D,S,T}(::Type{ScaledState{D,S,T}}) = D

	structure(::Type{ScaledState}) = AbstractStructure
	structure{D,S,T}(::Type{ScaledState{D,S,T}}) = S

	coeff(s::ScaledState) = s.coeff
	state(s::ScaledState) = s.state
	label(s::ScaledState) = label(state(s))

	######################
	# Printing Functions #
	######################
	repr(s::ScaledState) = "$(coeff(s)) $(repr(state(s)))"

	###########################
	# Mathematical Operations #
	###########################
	ctranspose(s::ScaledState) = ScaledState(coeff(s)', state(s)')

####################
# Both State Types #
####################
	###########################
	# Mathematical Operations #
	###########################
	# Note: inner is defined in 
	# scalar.jl, while outer is 
	# defined in diracoperator.jl
	
	kron{D,A,B}(a::DiracState{D, A}, b::DiracState{D, B}) = DiracState{D, typejoin(A,B)}(combine(label(a), label(b))) 
	kron{D}(a::AbstractDiracState{D}, b::AbstractDiracState{D}) = ScaledState(coeff(a)*coeff(b), kron(state(a), state(b))) 
	kron(a::AbstractDiracState, b::AbstractDiracState) = outer(a,b)

	for op=(:*, :.*)
		@eval begin
			($op){D}(a::AbstractDiracState{D}, b::AbstractDiracState{D}) = kron(a,b)
			($op)(c::Number, s::AbstractDiracState) = ScaledState(($op)(c, coeff(s)), state(s))
			($op)(s::AbstractDiracState, c::Number) = ScaledState(($op)(coeff(s), c), state(s))
			($op)(a::AbstractDiracBra, b::AbstractDiracKet) = inner(a, b)
			($op)(a::AbstractDiracKet, b::AbstractDiracBra) = outer(a, b)
		end
	end

	Ac_mul_B(a::AbstractDiracKet, b::AbstractDiracKet) = inner(a', b)
	Ac_mul_B(a::AbstractDiracBra, b::AbstractDiracBra) = outer(a', b)

	-(s::AbstractDiracState) = ScaledState(-(coeff(s)), state(s)) 
	+(s::AbstractDiracState) = s

export DualType,
	Ket,
	Bra,
	AbstractDiracState,
	AbstractDiracKet,
	AbstractDiracBra,
	DiracState,
	DiracKet,
	DiracBra,
	ScaledState,
	ScaledKet,
	ScaledBra,
	ket,
	bra,
	dualtype,
	structure,
	label,
	coeff,
	coefftype,
	state,
	samelabels
