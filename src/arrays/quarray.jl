###########
# QuArray #
###########
    abstract AbstractQuArray{B<:AbstractBasis,T,N}

    typealias AbstractQuVector{B<:AbstractBasis,T} AbstractQuArray{B,T,1}
    typealias AbstractQuMatrix{B<:AbstractBasis,T} AbstractQuArray{B,T,2}

    type QuArray{B<:AbstractBasis,T,N,A} <: AbstractQuArray{B,T,N}
        coeffs::A
        bases::NTuple{N,B}
        function QuArray(coeffs::AbstractArray{T,N}, bases::NTuple{N,B}) 
            if checkbases(coeffs, bases) 
                new(coeffs, bases)
            else 
                error("Coefficient array does not conform to input bases")
            end
        end
    end

    QuArray{B<:AbstractBasis,T,N}(coeffs::AbstractArray{T,N}, bases::NTuple{N,B}) = QuArray{B,T,N,typeof(coeffs)}(coeffs, bases)    
    QuArray(coeffs, bases::AbstractBasis...) = QuArray(coeffs, bases)
    QuArray(coeffs) = QuArray(coeffs, basesfordims(size(coeffs)))

    typealias QuVector{B<:AbstractBasis,T,N,A} QuArray{B,T,1,A}
    typealias QuMatrix{B<:AbstractBasis,T,N,A} QuArray{B,T,2,A}
 
    ######################
    # Accessor functions #
    ######################
    coefftype{B,T,N,A}(::QuArray{B,T,N,A}) = A
    coeffs(qarr::QuArray) = qarr.coeffs

    bases(qarr::QuArray, i) = qarr.bases[i]
    bases(qarr::AbstractQuArray) = ntuple(ndims(qarr), i->bases(qarr, i))

    ########################
    # Array-like functions #
    ########################
    Base.size(qarr::QuArray, i...) = size(coeffs(qarr), i...)

    Base.ndims(qarr::QuArray) = ndims(coeffs(qarr))
    Base.length(qarr::QuArray) = length(coeffs(qarr))

    Base.getindex(qarr::QuArray, i...) = getindex(coeffs(qarr), i...)
    Base.setindex!(qarr::QuArray, i...) = setindex!(coeffs(qarr), i...)

    Base.in(c, qarr::QuArray) = in(c, coeffs(qarr))

    ######################
    # Printing Functions #
    ######################
    Base.summary{B}(qarr::AbstractQuArray{B}) = "$(sizenotation(size(qarr))) QuArray in basis $B"

    function Base.show(io::IO, qarr::AbstractQuArray)
        println(io, summary(qarr)*":")
        println(io, "...transposed/conjugated?: $(istran(qarr))/$(isconj(qarr))")
        println(io, "...original coefficients: $(coefftype(qarr))")
        print(io, repr(coeffs(qarr)))
    end

    ####################
    # Helper Functions #
    ####################
    sizenotation(tup::(Int,)) = "$(first(tup))-element"
    sizenotation(tup::(Int...)) = reduce(*, map(s->"$(s)x", tup))[1:end-1] 

    function checkbases{N}(coeffs, bases::NTuple{N,AbstractBasis})
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

export QuArray
