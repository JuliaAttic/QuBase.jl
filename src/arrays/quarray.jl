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
    getbasis(qa::QuArray,i) = qa.bases[i]
    coefftype{B,T,N,A}(::QuArray{B,T,N,A}) = A
    coeffs(qa::QuArray) = qa.coeffs

    ########################
    # Array-like functions #
    ########################
    Base.size(qa::AbstractQuArray, i...) = size(coeffs(qa), i...)

    Base.ndims(qa::AbstractQuArray) = ndims(coeffs(qa))
    Base.length(qa::AbstractQuArray) = length(coeffs(qa))

    Base.getindex(qa::AbstractQuArray, i...) = getindex(coeffs(qa), i...)
    Base.setindex!(qa::AbstractQuArray, i...) = setindex!(coeffs(qa), i...)

    Base.in(c, qa::AbstractQuArray) = in(c, coeffs(qa))

    ######################
    # Printing Functions #
    ######################
    function Base.summary{B,T,N}(qa::AbstractQuArray{B,T,N})
        return "$(sizenotation(size(qa))) QuArray:\n" *
               "...bases: $B,\n" * 
               "...coefficients: $(typeof(coeffs(qa)))\n" *
               "...conjugate: $(isconj(qa))\n" *   
               "...transpose: $(istran(qa))"
    end

    # Right now just show the summary; we need to implement
    # array printing soon.
    Base.show(io::IO, qa::AbstractQuArray) = print(io, summary(qa))

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
