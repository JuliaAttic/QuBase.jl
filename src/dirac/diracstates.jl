import Base: getindex,
    length,
    show,
    repr,
    convert,
    ctranspose,
    Ac_mul_B, #TODO: Define other members of the A_mul_B family as necessary 
    +,
    -,
    *,.*

################################
# Abstract types and Functions #
################################
    # An AbstractState is a type that
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
    # Concrete subtypes of AbstractState 
    # implement the following methods (T<:AbstractState):
    #
    # coeff(s::T) -> returns the scalar coefficient of the state
    # state(s::T) -> returns only the state, without a scalar coefficient, as a DiracState
    # label(s::T) -> returns the state's label as a StateLabel
    # repr(s::T) -> returns a string representation of the state
    # coefftype(s::T), 
    # coefftype(::Type{T}) -> returns the coefficient type for a given state/type

    getindex(s::AbstractState, i) = getindex(label(s), i)
    length(s::AbstractState) = length(label(s))
    show(io::IO, s::AbstractState) = print(io, repr(s))

    # This comparison method is useful to define
    # for all labeled concrete subtypes of AbstractQuantum
    samelabels(a::AbstractState, b::AbstractState) = label(a) == label(b)

    convert(::Type{StateLabel}, s::AbstractState) = label(s)

##############
# DiracState #
##############
    # A DiracState is the type representation of
    # an unscaled abstract vector, formulated in 
    # Dirac notation as a bra or ket.
    immutable DiracState{D<:DualType, S<:AbstractStructure} <: AbstractState{D, S}
        label::StateLabel
    end

    typealias DiracKet{S} DiracState{Ket, S}
    typealias DiracBra{S} DiracState{Bra, S}

    convert{D,S}(::Type{DiracState{D,S}}, s::AbstractState) = DiracState{D,S}(label(s))

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
    # AbstractState Functions #
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
    ctranspose{D,S}(k::DiracState{D,S}) = DiracState{D',S}(label(k))

###############
# ScaledState #
###############
    # A ScaledState is the type representation of
    # a DiracState multiplied by a scalar. 
    immutable ScaledState{D<:DualType, S<:AbstractStructure, T} <: AbstractState{D, S}
        coeff::T
        state::DiracState{D,S}
    end
        
    convert{D,S,T}(::Type{ScaledState{D,S,T}}, s::AbstractState) = ScaledState(convert(T, coeff(s)), convert(DiracState{D,S}, state(s)))

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
    
    tensor{D,A,B}(a::DiracState{D, A}, b::DiracState{D, B}) = DiracState{D, typejoin(A,B)}(combine(label(a), label(b))) 
    tensor{D}(a::AbstractState{D}, b::AbstractState{D}) = ScaledState(kron(coeff(a),coeff(b)), tensor(state(a), state(b))) 
    tensor(a::AbstractState{Ket}, b::AbstractState{Bra}) = outer(a,b)
    tensor(a::AbstractState{Bra}, b::AbstractState{Ket}) = outer(b,a)

    for op=(:*, :.*)
        @eval begin
            ($op){D}(a::AbstractState{D}, b::AbstractState{D}) = tensor(a,b)
            ($op)(c::Number, s::AbstractState) = ScaledState(($op)(c, coeff(s)), state(s))
            ($op)(s::AbstractState, c::Number) = ScaledState(($op)(coeff(s), c), state(s))
            ($op)(a::AbstractBra, b::AbstractKet) = inner(a, b)
            ($op)(a::AbstractKet, b::AbstractBra) = outer(a, b)
        end
    end

    Ac_mul_B(a::AbstractKet, b::AbstractKet) = inner(a', b)
    Ac_mul_B(a::AbstractBra, b::AbstractBra) = outer(a', b)

    -(s::AbstractState) = ScaledState(-(coeff(s)), state(s)) 
    +(s::AbstractState) = s

export DiracState,
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
    samelabels,
    tensor
