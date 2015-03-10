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

    immutable FiniteBasis{S,N} <: AbstractFiniteBasis{S}
        lens::NTuple{N,Int}
        FiniteBasis(lens::NTuple{N,Int}) = new(lens)
        FiniteBasis(lens::Int...) = new(lens)
    end

    FiniteBasis{N, S<:AbstractStructure}(lens::NTuple{N,Int}, ::Type{S}=Orthonormal) = FiniteBasis{S,N}(lens)
    FiniteBasis(lens::Int...) = FiniteBasis(lens)

    Base.convert{S,N}(::Type{FiniteBasis{S,N}}, f::FiniteBasis) = FiniteBasis(f.lens, S)

    ######################
    # Property Functions #
    ######################
    structure{S}(::Type{FiniteBasis{S}}) = S
    structure{S,N}(::Type{FiniteBasis{S,N}}) = S

    Base.length(basis::FiniteBasis) = prod(basis.lens)
    Base.size(basis::FiniteBasis) = basis.lens
    Base.size(basis::FiniteBasis, i) = basis.lens[i]
    nfactors{S,N}(::FiniteBasis{S,N}) = N
    checkcoeffs(coeffs, dim::Int, basis::FiniteBasis) = size(coeffs, dim) == length(basis) 

    ##########################
    # Mathematical Functions #
    ##########################
    tensor{S}(a::FiniteBasis{S}, b::FiniteBasis{S}) = FiniteBasis(tuple(a.lens..., b.lens...), S)

    ######################
    # Printing Functions #
    ######################
    Base.repr(f::FiniteBasis) = "$(typeof(f))$(size(f))"
    Base.show(io::IO, f::FiniteBasis) = print(io, repr(f))

export FiniteBasis,
    structure,
    tensor,
    checkcoeffs