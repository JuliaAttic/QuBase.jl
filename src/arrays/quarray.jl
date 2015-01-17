import Base: size,
    length,
    getindex,
    similar,
    in,
    ctranspose,
    transpose,
    summary,
    #TODO: Implement the below operations
    *,.*,
    /,./,
    +,.+,
    -,.-,
    kron

abstract AbstractQuArray{B<:AbstractBasis,T,N} <: AbstractArray{T,N}

typealias AbstractQuVector{B<:AbstractBasis,T} AbstractQuArray{B,T,1}
typealias AbstractQuMatrix{B<:AbstractBasis,T} AbstractQuArray{B,T,2}

###########
# QuArray #                 
###########
    # Using an NTuple allows us to have a basis for each dimension, 
    # and gives us a less ambiguous way to determine if a QuArray will act
    # like a vector, matrix, or tensor using the dimension parameter N. 
    type QuArray{B<:AbstractBasis,T,N,A} <: AbstractQuArray{B,T,N}
        coeffs::A
        bases::NTuple{N,B}
        function QuArray(coeffs::AbstractArray{T}, bases::NTuple{N,B}) 
            if checkbases(coeffs, bases) 
                new(coeffs, bases)
            else 
                error("Coefficient array does not conform to input bases")
            end
        end
    end
    
    typealias QuVector{B<:AbstractBasis,T,A} QuArray{B,T,1,A}
    typealias QuMatrix{B<:AbstractBasis,T,A} QuArray{B,T,2,A}

    QuArray{T,N,B<:AbstractBasis}(coeffs::AbstractArray{T}, bases::NTuple{N,B}) = QuArray{B,T,N,typeof(coeffs)}(coeffs, bases)
    QuArray(coeffs, bases::AbstractBasis...) = QuArray(coeffs, bases)
    QuArray(coeffs) = QuArray(coeffs, basesfordims(size(coeffs)))
    
    ######################
    # Property Functions #
    ######################
    bases(qa::QuArray) = qa.bases
    coeffs(qa::QuArray) = qa.coeffs
    size(qa::QuArray, i...) = size(qa.coeffs, i...)

    ########################
    # Array-like functions #
    ########################
    similar{B,T}(qa::QuArray{B,T}, element_type=T) = QuArray(similar(qa.coeffs, T), qa.bases)
    # Is there a way to properly define the below for
    # any arbitrary basis? Obviously doesn't make sense
    # for B<:AbstractInfiniteBasis, and I'm reluctant to
    # enforce that every B<:AbstractFiniteBasis will have a 
    # constructor B(::Int), which is how the below is constructing
    # instances of FiniteBasis.
    function similar{B<:FiniteBasis,T}(qa::QuArray{B,T}, element_type, dims)
        return QuArray(similar(qa.coeffs, T, dims), basesfordims(dims, B))
    end 
    getindex(qa::QuArray, i::AbstractArray) = getindex(qa.coeffs, i)
    getindex(qa::QuArray, i::Real) = getindex(qa.coeffs, i)
    getindex(qa::QuArray, i) = getindex(qa.coeffs, i)
    getindex(qa::QuArray, i...) = getindex(qa.coeffs, i...)

    in(c, qa::QuArray) = in(c, qa.coeffs)

    ctranspose(qa::QuVector) = QuArray(ctranspose(qa.coeffs), qa.bases)
    ctranspose(qa::QuMatrix) = QuArray(ctranspose(qa.coeffs), reverse(qa.bases))
    transpose(qa::QuVector) = QuArray(transpose(qa.coeffs), qa.bases)
    transpose(qa::QuMatrix) = QuArray(transpose(qa.coeffs), reverse(qa.bases))

    ######################
    # Printing Functions #
    ######################
    function summary{B,T,N,A}(qa::QuArray{B,T,N,A})
        return "$(sizenotation(size(qa))) QuArray\n" *
               "...bases: $B,\n" * 
               "...coeff: $A"             
    end

    ######################
    # Include Statements #
    ######################
    include("constructors.jl")
    include("ladderops.jl")
    
    ####################
    # Helper Functions #
    ####################
    sizenotation(tup::(Int,)) = "$(first(tup))-element"
    sizenotation(tup::(Int...)) = reduce(*, map(s->"$(s)x", tup))[1:end-1] 

    # checkbases() is overloaded for a single basis
    # to handle the fact that row vectors are
    # N=2
    function checkbases(coeffs, bases::NTuple{1, AbstractBasis})
        if ndims(coeffs) <= 2
            if size(coeffs, 2) == 1
                return checkcoeffs(coeffs, 1, first(bases))
            elseif size(coeffs, 1) == 1
                return checkcoeffs(coeffs, 2, first(bases))
            end 
        end
        return false
    end

    function checkbases{N}(coeffs, bases::NTuple{N, AbstractBasis})
        if ndims(coeffs) == length(bases)
            return reduce(&, [checkcoeffs(coeffs, i, bases[i]) for i=1:N])
        end
        return false
    end

    # Assumes that every basis type passed in
    # has a constructor B(::eltype(lens))
    function basesfordims(lens::Tuple, B=ntuple(length(lens), x->FiniteBasis))
        return ntuple(length(lens), n->B[n](lens[n]))
    end
    one_at_ind!(arr, i) = setindex!(arr, one(eltype(arr)), i)
    single_coeff(i, lens...) = one_at_ind!(zeros(lens), i)

export AbstractQuArray,
    QuArray,
    QuVector,
    QuMatrix,
    bases,
    coeffs
