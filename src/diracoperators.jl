import Base:
	show,
	repr,
	convert,
	ctranspose,
	Ac_mul_B, #TODO: Define other members of the A_mul_B family as necessary 
	+,
	-

################################
# Abstract types and Functions #
################################
	# An AbstractDiracOperator is a type that
	# represents an operator formulated in
	# Dirac notation - a ket-bra projector 
	# optionally multiplied by a scalar value. 
	#
	# Concrete subtypes of AbstractDiracOperator
	# implement the following methods:
	#
	# coeff(op::AbstractDiracOperator) -> returns the scalar coefficient of the operator
	# operator(op::AbstractDiracOperator) -> returns only the bra-ket projector, without a scalar coefficient, as a DiracOperator
	# getbra(op::AbstractDiracOperator) -> returns the bra component of the operator as a DiracBra
	# getket(op::AbstractDiracOperator) -> returns the ket component of the operator as a DiracKet
	# label(op::AbstractDiracOperator) -> returns the operator's labels as a tuple of type (StateLabel, StateLabel)
	# repr(op::AbstractDiracOperator) -> returns a string representation of the operator
	# coefftype(op::AbstractDiracOperator), 
	# coefftype(::Type{AbstractDiracOperator}) -> returns the coefficient type for a given operator/type
	# structure(::Type{AbstractDiracOperator}) -> returns the operator type's structure type parameter

	abstract AbstractDiracOperator{S<:AbstractStructure} <: AbstractOperator{S}

	structure{S}(::Type{AbstractDiracOperator{S}}) = S
	structure(::Type{AbstractDiracOperator}) = AbstractStructure

	show(io::IO, op::AbstractDiracOperator) = print(io, repr(op))

	samelabels(a::AbstractDiracOperator, b::AbstractDiracOperator) = samelabels(getket(a), getket(b)) && samelabels(getbra(a), getbra(b))

#################
# DiracOperator #
#################
	# A DiracOperator is the type representation of
	# an unscaled abstract operator, formulated in 
	# Dirac notation as a ket-bra projector.
	immutable DiracOperator{S} <: AbstractDiracOperator{S}
		ket::DiracKet{S}
		bra::DiracBra{S}
	end

	DiracOperator{A,B}(k::DiracKet{A}, b::DiracBra{B}) = DiracOperator{typejoin(A,B)}(k, b)

	getket(op::DiracOperator) = op.ket
	getbra(op::DiracOperator) = op.bra

	######################
	# Accessor Functions #
	######################
	coefftype(op::DiracOperator) = Int
	coefftype{O<:DiracOperator}(::Type{O}) = Int

	structure(::Type{DiracOperator}) = AbstractStructure
	structure{S}(::Type{DiracOperator{S}}) = S

	coeff(op::DiracOperator) = one(coefftype(op))
	operator(op::DiracOperator) = op
	label(op::DiracOperator) = (label(getket(op)), label(getbra(op)))

	######################
	# Printing Functions #
	######################
	repr(op::DiracOperator) = repr(getket(op))*repr(getbra(op))

	###########################
	# Mathematical Operations #
	###########################
	outer(k::DiracKet, b::DiracBra) = DiracOperator(k,b)
	ctranspose(op::DiracOperator) = outer(getbra(op)', getket(op)')

##################
# ScaledOperator #
##################
	# A ScaledOperator is the type representation of
	# a DiracOperator multiplied by a scalar.
	immutable ScaledOperator{S<:AbstractStructure, T} <: AbstractDiracOperator{S}
		coeff::T
		operator::DiracOperator{S}
	end

	######################
	# Accessor Functions #
	######################
	coefftype{S,T}(::ScaledOperator{S,T}) = T
	coefftype(::Type{ScaledOperator}) = Any
	coefftype{S,T}(::Type{ScaledOperator{S,T}}) = T

	structure(::Type{ScaledOperator}) = AbstractStructure
	structure{S,T}(::Type{ScaledOperator{S,T}}) = S

	coeff(op::ScaledOperator) = op.coeff
	operator(op::ScaledOperator) = op.operator
	label(op::ScaledOperator) = label(operator(op))

	######################
	# Printing Functions #
	######################
	repr(op::ScaledOperator) = "$(coeff(op)) $(repr(operator(op)))"

	###########################
	# Mathematical Operations #
	###########################
	ctranspose(op::ScaledOperator) = ScaledOperator(coeff(op)', operator(op)')
	outer(k::AbstractDiracKet, b::AbstractDiracBra) = ScaledOperator(kron(coeff(k), coeff(b)), outer(state(k), state(b)))

#######################
# Both Operator Types #
#######################
	###########################
	# Mathematical Operations #
	###########################
	tensor(a::DiracOperator, b::DiracOperator) = outer(tensor(getket(a), getket(b)), tensor(getbra(a), getbra(b)))
	tensor(a::AbstractDiracOperator, b::AbstractDiracOperator) = ScaledOperator(kron(coeff(a),coeff(b)), tensor(operator(a), operator(b))) 

	tensor(op::DiracOperator, s::DiracBra) = outer(getket(op), tensor(getbra(op), s))
	tensor(s::DiracBra, op::DiracOperator) = outer(getket(op), tensor(s, getbra(op)))
	tensor(op::DiracOperator, s::DiracKet) = outer(tensor(getket(op), s), getbra(op))
	tensor(s::DiracKet, op::DiracOperator) = outer(tensor(s, getket(op)), getbra(op))
	tensor(op::AbstractDiracOperator, s::AbstractDiracState) = ScaledOperator(kron(coeff(op),coeff(s)), tensor(operator(op), s)) 
	tensor(s::AbstractDiracState, op::AbstractDiracOperator) = ScaledOperator(kron(coeff(s),coeff(op)), tensor(s, operator(op)))

	for op=(:*, :.*)
		@eval begin
			($op)(s::AbstractDiracKet, op::AbstractDiracOperator) = tensor(s,op)
			($op)(op::AbstractDiracOperator, s::AbstractDiracBra) = tensor(op,s)
			($op)(a::AbstractDiracOperator, b::AbstractDiracOperator) = ScaledOperator(inner(getbra(a), getket(b)), tensor(getket(a), getbra(b))) 
			($op)(s::AbstractDiracBra, op::AbstractDiracOperator) = ($op)(inner(s, getket(op)), getbra(op))
			($op)(op::AbstractDiracOperator, s::AbstractDiracKet) = ($op)(inner(getket(op), s), getket(op))
			($op)(c::Number, op::AbstractDiracOperator) = ScaledOperator(($op)(c,coeff(op)), operator(op))
			($op)(op::AbstractDiracOperator, c::Number) = ScaledOperator(($op)(coeff(op),c), operator(op))
		end
	end

	Ac_mul_B(a::AbstractDiracOperator, b::AbstractDiracOperator) = inner(a', b)
	Ac_mul_B(a::AbstractDiracBra, b::AbstractDiracOperator) = tensor(a', b)
	Ac_mul_B(a::AbstractDiracOperator, b::AbstractDiracBra) = tensor(a', b)
	Ac_mul_B(a::AbstractDiracKet, b::AbstractDiracOperator) = inner(a', getket(b)) * getbra(a)
	Ac_mul_B(a::AbstractDiracOperator, b::AbstractDiracKet) = inner(getket(a)', b) * getbra(a)'

	-(op::AbstractDiracOperator) = ScaledOperator(-(coeff(op)), operator(op)) 
	+(op::AbstractDiracOperator) = op

export AbstractDiracOperator,
	DiracOperator,
	ScaledOperator,
	outer,
	operator,
	getket,
	getbra,
	coeff,
	coefftype,
	label,
	samelabels,
	structure,
	tensor