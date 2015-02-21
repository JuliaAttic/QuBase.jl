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
            if N <= 2
                if checkbases(coeffs, bases) 
                    new(coeffs, bases)
                else 
                    error("Coefficient array does not conform to input bases")
                end
            else
                error("QuArrays are only defined for coefficient arrays of dimension N<=2")
            end
        end
    end

    QuArray{B<:AbstractBasis,T,N}(coeffs::AbstractArray{T,N}, bases::NTuple{N,B}) = QuArray{B,T,N,typeof(coeffs)}(coeffs, bases)    
    QuArray(coeffs, bases::AbstractBasis...) = QuArray(coeffs, bases)
    QuArray(coeffs) = QuArray(coeffs, basesfordims(size(coeffs)))
    
    #####################################
    # Conjugate/Transpose Wrapper types #
    #####################################
    immutable Transpose{B,T,N,A,Q<:AbstractQuArray{B,T,N}} <: AbstractQuArray{B,T,N}
        qa::Q    
    end

    immutable Conjugate{B,T,N,A,Q<:AbstractQuArray{B,T,N}} <: AbstractQuArray{B,T,N}
        qa::Q
    end

    Transpose{B,T,N}(qa::AbstractQuArray{B,T,N}) = Transpose{B,T,N,typeof(coeffs(qa)),typeof(qa)}(qa)
    Conjugate{B,T,N}(qa::AbstractQuArray{B,T,N}) = Conjugate{B,T,N,typeof(coeffs(qa)),typeof(qa)}(qa)

    #######################
    # QuArray typealiases #
    #######################
    typealias ConjTran{B<:AbstractBasis,T,N,A,Q<:Transpose} Conjugate{B,T,N,A,Q}
    typealias TranConj{B<:AbstractBasis,T,N,A,Q<:Conjugate} Transpose{B,T,N,A,Q}
    
    typealias TranArray{B<:AbstractBasis,T,N,A} Union(ConjTran{B,T,N,A}, Transpose{B,T,N,A}) 
    typealias TranVector{B<:AbstractBasis,T,N,A} TranArray{B,T,1,A}
    typealias TranMatrix{B<:AbstractBasis,T,N,A} TranArray{B,T,2,A}

    typealias ConjArray{B<:AbstractBasis,T,N,A} Union(TranConj{B,T,N,A}, Conjugate{B,T,N,A}) 
    typealias ConjVector{B<:AbstractBasis,T,N,A} ConjArray{B,T,1,A}
    typealias ConjMatrix{B<:AbstractBasis,T,N,A} ConjArray{B,T,2,A}

    typealias DualArray{B<:AbstractBasis,T,N,A} Union(ConjTran{B,T,N,A}, TranConj{B,T,N,A})
    typealias DualVector{B<:AbstractBasis,T,N,A} DualArray{B,T,1,A}
    typealias DualMatrix{B<:AbstractBasis,T,N,A} DualArray{B,T,2,A}

    typealias QuVector{B<:AbstractBasis,T,N,A} QuArray{B,T,1,A}
    typealias QuMatrix{B<:AbstractBasis,T,N,A} QuArray{B,T,2,A}

    ######################
    # Accessor functions #
    ######################
    qarr(qa::QuArray) = qa
    qarr(qt::Transpose) = qarr(qt.qa)
    qarr(qc::Conjugate) = qarr(qc.qa)
    coeffs(qa::AbstractQuArray) = qarr(qa).coeffs
    
    ########################
    # Array-like functions #
    ########################
    Base.size(qa::AbstractQuArray, i...) = size(coeffs(qa), i...)
    Base.size(qt::TranMatrix) = reverse(size(coeffs(qt)))
    Base.size(qt::TranMatrix, i) = size(qt)[i]

    Base.ndims(qa::AbstractQuArray) = ndims(coeffs(qa))
    Base.length(qa::AbstractQuArray) = length(coeffs(qa))

    getbasis(qa::AbstractQuArray,i) = qarr(qa).bases[i]
    getbasis(qt::TranMatrix,i) = qarr(qa).bases[(length+1)-i]

    Base.getindex(qa::QuArray, i...) = getindex(qa.coeffs, i...)
    Base.getindex(qc::Conjugate, i...) = conj(qc.qa[i...])
    Base.getindex(qt::Transpose, x, i, j, k...) = qt.qa[j,i,k...]
    Base.getindex(qt::TranVector, i) = qt.qa[i]

    Base.setindex!(qa::QuArray, i) = setindex!(qa.coeffs, i)
    Base.setindex!(qa::QuArray, i...) = setindex!(qa.coeffs, i...)
    Base.setindex!(qc::Conjugate, x, i...) = (qc.qa[i...] = conj(x))
    Base.setindex!(qt::Transpose, x, i, j, k...) = (qt.qa[j,i,k...] = x)
    Base.setindex!(qt::TranMatrix, x, i) = (qt.qa[i] = x)

    Base.in(c, qa::AbstractQuArray) = in(c, coeffs(qa))

    Base.conj(qa::AbstractQuArray) = Conjugate(qa)
    Base.conj(qt::Transpose) = Transpose(conj(qt.qa))
    Base.conj(qc::Conjugate) = qc.qa

    Base.transpose(qa::AbstractQuArray) = Transpose(qa)
    Base.transpose(qc::Conjugate) = Conjugate(transpose(qc.qa))
    Base.transpose(qt::Transpose) = qt.qa
    Base.ctranspose(qa::AbstractQuArray) = conj(transpose(qa))

    ######################
    # Printing Functions #
    ######################
    function Base.summary{B,T,N}(qa::AbstractQuArray{B,T,N})
        return "$(sizenotation(size(qa))) QuArray:\n" *
               "...bases: $B,\n" * 
               "...coefficients: $(typeof(coeffs(qa)))\n" *
               "...conjugate: $(typeof(qa) <: ConjArray)\n" *   
               "...transpose: $(typeof(qa) <: TranArray)"
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
