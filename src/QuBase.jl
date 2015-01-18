module QuBase
    
    import Base: kron
    
    ####################
    # String Constants #
    ####################
        const lang = "\u27E8"
        const rang = "\u27E9"
        const otimes = "\u2297"
        const vdots ="\u205E"

    ##################
    # Abstract Types #
    ##################
        abstract AbstractStructure
        abstract Orthonormal <: AbstractStructure
        abstract AbstractQuantum{S<:AbstractStructure}

        # Various constructor methods in this repo allow an argument 
        # of type Type{BypassFlag} to be passed in in order to 
        # circumvent value precalculation/checking. This is useful for
        # conversion methods and the like, where you know the input 
        # has already been vetted elsewhere. Don't use this unless
        # you're sure of what you're doing, and don't export this.
        abstract BypassFlag

    #############
    # Functions #
    #############

        # This should be the only `structure` 
        # method that needs to be defined for 
        # type *instances*
        structure{S}(::AbstractQuantum{S}) = S
        structure(::DataType) = AbstractStructure

        # ...and all relevant singleton types 
        # should have it defined as well:
        structure{S}(::Type{AbstractQuantum{S}}) = S

        # an n-arity form of the tensor
        # product, reduction is done via
        # the binary definition of tensor()
        # defined in the files included above. 
        tensor(s...) = reduce(tensor, s)

        # For the sake of convenience, kron() 
        # is defined to be equivalent to 
        # tensor() for quantum objects
        kron(a::AbstractQuantum, b::AbstractQuantum) = tensor(a, b)
    
    ######################
    # Include Statements #
    ######################
        include("bases/bases.jl")
        include("arrays/arrays.jl")
    
    export AbstractStructure, 
        AbstractQuantum,
        Orthonormal,
        structure,
        tensor 
end
    