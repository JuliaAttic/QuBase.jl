###########
# QuArray #
###########
    abstract AbstractQuArray{B<:AbstractBasis,T,N}

    typealias AbstractQuVector{B<:AbstractBasis,T} AbstractQuArray{B,T,1}
    typealias AbstractQuMatrix{B<:AbstractBasis,T} AbstractQuArray{B,T,2}

    type QuArray{B<:AbstractBasis,T,N,C} <: AbstractQuArray{B,T,N}
        coeffs::C
        bases::NTuple{N,B}
        function QuArray{Tran,Conj}(coeffs::QuCoeffs{Tran,Conj,N,T}, bases::NTuple{N,B}) 
            if checkbases(coeffs, bases) 
                new(coeffs, bases)
            else 
                error("Coefficient array does not conform to input bases")
            end
        end
    end
    
    typealias QuVector{B<:AbstractBasis,T,C} QuArray{B,T,1,C}
    typealias QuMatrix{B<:AbstractBasis,T,C} QuArray{B,T,2,C}

    typealias QuKet{B<:AbstractBasis,T,KC<:KetCoeffs} QuVector{B,T,KC}
    typealias QuBra{B<:AbstractBasis,T,BC<:BraCoeffs} QuVector{B,T,BC}

    function QuArray{Tran,Conj,T,N,B<:AbstractBasis}(coeffs::QuCoeffs{Tran,Conj,N,T}, 
                                                     bases::NTuple{N,B})
        return QuArray{B,T,N,typeof(coeffs)}(coeffs, bases)
    end

    QuArray{N,B<:AbstractBasis}(coeffs::AbstractArray, bases::NTuple{N,B}) = QuArray(QuCoeffs(coeffs), bases)
    QuArray(coeffs, bases::AbstractBasis...) = QuArray(coeffs, bases)
    QuArray(coeffs) = QuArray(coeffs, basesfordims(size(coeffs)))
    
    ######################
    # Property Functions #
    ######################
    bases(qa::QuArray) = qa.bases
    coeffs(qa::QuArray) = qa.coeffs

    ########################
    # Array-like functions #
    ########################
    Base.size(qa::AbstractQuArray, i...) = size(coeffs(qa), i...)
    Base.ndims(qa::AbstractQuArray) = ndims(coeffs(qa))
    Base.length(qa::AbstractQuArray) = length(coeffs(qa))

    Base.getindex(qa::AbstractQuArray, i) = getindex(coeffs(qa), i)
    Base.getindex(qa::AbstractQuArray, i...) = getindex(coeffs(qa), i...)

    Base.setindex!(qa::AbstractQuArray, i) = setindex!(coeffs(qa), i)
    Base.setindex!(qa::AbstractQuArray, i...) = setindex!(coeffs(qa), i...)

    Base.in(c, qa::AbstractQuArray) = in(c, coeffs(qa))

    Base.conj(qa::AbstractQuArray) = QuArray(conj(coeffs(qa)), bases(qa))
    Base.transpose(qa::AbstractQuArray) = QuArray(transpose(coeffs(qa)), reverse(bases(qa)))
    Base.ctranspose(qa::AbstractQuArray) = QuArray(ctranspose(coeffs(qa)), reverse(bases(qa)))

    ######################
    # Printing Functions #
    ######################
    function Base.summary{B,T,N,C}(qa::QuArray{B,T,N,C})
        return "$(sizenotation(size(qa))) QuArray\n" *
               "...bases: $B,\n" * 
               "...coeff: $C"             
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

export AbstractQuArray,
    QuArray,
    QuVector,
    QuMatrix,
    bases,
    coeffs
