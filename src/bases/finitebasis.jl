import Base: length,
    size,
    ndims

###############
# FiniteBasis #
###############
    # A FiniteBasis repesents the tensor product
    # of an arbitrary number of fock bases, each of which is truncated 
    # at a known, finite number of states. For example, imagine the 
    # following fock bases and their tensor product:
    #
    # F_4 => {|0>, |1>, |2>, |3>}
    # F_3 => {|0>, |1>, |2>}
    # F_2 => {|0>, |1>} 
    #
    # F_4 x F_3 x F_2 => {|0,0,0>, |0,0,1>, |0,1,0>, |0,1,1>...|3,2,1>}
    # 
    # Note that the subscript number refers to basis *length*, not max
    # excitation level.
    #
    # This tensor product basis, F_4 x F_3 x F_2, could be represented 
    # by calling FiniteBasis(4, 3, 2). 

    immutable FiniteBasis{S} <: AbstractBasis{S}
        lens::(Int...)
        FiniteBasis(lens::(Int...)) = new(lens)
        FiniteBasis(lens::Int...) = new(lens)
    end

    FiniteBasis(lens::(Int...)) = FiniteBasis{AbstractStructure}(lens)
    FiniteBasis(lens::Int...) = FiniteBasis{AbstractStructure}(lens)

    convert{S}(::Type{FiniteBasis{S}}, f::FiniteBasis) = FiniteBasis{S}(f.lens)

    ######################
    # Property Functions #
    ######################
    structure{S}(::Type{FiniteBasis{S}}) = S
    structure(::Type{FiniteBasis}) = AbstractStructure

    length(basis::FiniteBasis) = prod(basis.lens)
    size(basis::FiniteBasis) = basis.lens
    size(basis::FiniteBasis, i) = basis.lens[i]
    ndims(basis::FiniteBasis) = length(basis.lens)

    checkcoeffs(coeffs::AbstractArray, dim::Int, basis::FiniteBasis) = size(coeffs, dim) == length(basis) 

    ##########################
    # Mathematical Functions #
    ##########################
    tensor{A,B}(a::FiniteBasis{A}, b::FiniteBasis{B}) = FiniteBasis{typejoin(A,B)}(tuple(a.lens..., b.lens...))

    ######################
    # Printing Functions #
    ######################
    repr(f::FiniteBasis) = "$(typeof(f))$(size(f))"
    show(io::IO, f::FiniteBasis) = print(io, repr(f))

export FiniteBasis,
    structure,
    checkcoeffs