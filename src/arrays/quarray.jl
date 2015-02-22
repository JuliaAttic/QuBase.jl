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

    QuArray{T,N,B<:AbstractBasis}(coeffs::AbstractArray{T,N}, bases::NTuple{N,B}) = QuArray{B,T,N,typeof(coeffs)}(coeffs, bases)
    QuArray(coeffs, bases::AbstractBasis...) = QuArray(coeffs, bases)
    QuArray(coeffs) = QuArray(coeffs, basesfordims(size(coeffs)))
    
    typealias QuCol{B<:AbstractBasis,T,A} QuArray{B,T,1,A}
    typealias QuMatrix{B<:AbstractBasis,T,A} QuArray{B,T,2,A}

    immutable QuRow{B,T,A} <: AbstractQuVector{B,T}
        qcol::QuCol{B,T,A}  
    end

    QuRow{B,T,A}(qcol::QuCol{B,T,A}) = QuRow{B,T,A}(qcol)

    ######################
    # Property Functions #
    ######################
    bases(qarr::QuArray) = qarr.bases
    coeffs(qarr::QuArray) = qarr.coeffs
    
    bases(qrow::QuRow) = bases(qrow.qcol)
    coeffs(qrow::QuRow) = coeffs(qrow.qcol)

    coeff_apply(f, qarr::QuArray) = QuArray(f(coeffs(qarr)), bases(qarr))
    coeff_apply(f, qrow::QuRow) = QuRow(coeff_apply(f, qrow.qcol))

    ########################
    # Array-like functions #
    ########################
    Base.size(qarr::AbstractQuArray, i...) = size(coeffs(qarr), i...)
    Base.ndims(qarr::AbstractQuArray) = ndims(coeffs(qarr))
    Base.length(qarr::AbstractQuArray) = length(coeffs(qarr))

    Base.getindex(qarr::AbstractQuArray, i) = getindex(coeffs(qarr), i)
    Base.getindex(qarr::AbstractQuArray, i...) = getindex(coeffs(qarr), i...)

    Base.setindex!(qarr::AbstractQuArray, i) = setindex!(coeffs(qarr), i)
    Base.setindex!(qarr::AbstractQuArray, i...) = setindex!(coeffs(qarr), i...)

    Base.in(c, qarr::AbstractQuArray) = in(c, coeffs(qarr))

    Base.conj(qarr::AbstractQuArray) = coeff_apply(conj, qarr)
    
    Base.transpose(qcol::QuCol) = QuRow(qcol)
    Base.transpose(qrow::QuRow) = qrow.qcol 
    Base.transpose(qarr::QuArray) = QuArray(transpose(coeffs(qarr)), reverse(bases(qarr)))
    
    Base.ctranspose(qcol::QuCol) = QuRow(conj(qcol))
    Base.ctranspose(qrow::QuRow) = conj(transpose(qrow))
    Base.ctranspose(qarr::QuArray) = QuArray(ctranspose(coeffs(qarr)), reverse(bases(qarr)))

    ######################
    # Printing Functions #
    ######################

    Base.summary{B}(qarr::AbstractQuArray{B}) = "$(sizenotation(size(qarr))) $(typenotation(qarr)) in basis $B"              

    function Base.show(io::IO, qarr::AbstractQuArray)
        println(io, summary(qarr)*":")
        println(io, "...coefficients: $(typeof(coeffs(qarr)))")
        print(io, repr(coeffs(qarr)))
    end

    ####################
    # Helper Functions #
    ####################
    typenotation(::QuArray) = "QuArray"
    typenotation(::QuMatrix) = "QuMatrix"
    typenotation(::QuRow) = "QuRow"
    typenotation(::QuCol) = "QuCol"

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
