############
# QuCoeffs #
############
    abstract ConjBool{Conj}
    abstract TranBool{Tran}

    # enumerate cases for type-stable flipping
    flip(::Type{ConjBool{true}}) = ConjBool{false}
    flip(::Type{TranBool{true}}) = TranBool{false}
    flip(::Type{ConjBool{false}}) = ConjBool{true}
    flip(::Type{TranBool{false}}) = TranBool{true}

    type QuCoeffs{Tran,Conj,N,T,A}
        arr::A  
        tran::Type{TranBool{Tran}}
        conj::Type{ConjBool{Conj}}
        function QuCoeffs(arr::AbstractArray{T}, 
                          tran::Type{TranBool{Tran}},
                          conj::Type{ConjBool{Conj}})
            return new(arr, tran, conj)
        end
    end

    function QuCoeffs{Conj,Tran,T,N}(arr::AbstractArray{T,N}, 
                      tran::Type{TranBool{Tran}},
                      conj::Type{ConjBool{Conj}})
        return QuCoeffs{Tran,Conj,N,T,typeof(arr)}(arr, tran, conj)
    end

    QuCoeffs(arr::AbstractArray) = QuCoeffs(arr, TranBool{false}, ConjBool{false})

    typealias StateCoeffs{Tran,Conj,T,A} QuCoeffs{Tran,Conj,1,T,A}
    typealias OpCoeffs{Tran,Conj,T,A} QuCoeffs{Tran,Conj,2,T,A}

    typealias KetCoeffs{T,A} StateCoeffs{false,false,T,A}
    typealias BraCoeffs{T,A} StateCoeffs{true,true,T,A}

    ########################
    # Array-like Functions #
    ########################
    Base.size(qc::QuCoeffs) = size(qc.arr)
    Base.size(qc::QuCoeffs, i) = size(qc.arr, i)

    Base.ndims(qc::QuCoeffs) = ndims(qc.arr)    
    Base.length(qc::QuCoeffs) = length(qc.arr)

    Base.getindex(qc::QuCoeffs, i...) = getindex(qc.arr, i...)
    Base.setindex!(qc::QuCoeffs, i...) = setindex!(qc.arr, i...)

    #######################
    # Conjugate/Transpose #
    #######################
    Base.conj(qc::QuCoeffs) = QuCoeffs(conj(qc.arr), qc.tran, flip(qc.conj))
    
    Base.transpose(qc::QuCoeffs) = QuCoeffs(transpose(qc.arr), flip(qc.tran), qc.conj)
    Base.transpose(qc::StateCoeffs) = QuCoeffs(copy(qc.arr), flip(qc.tran), qc.conj)
 
    Base.ctranspose(qc::QuCoeffs) = QuCoeffs(ctranspose(qc.arr), flip(qc.tran), flip(qc.conj))
    Base.ctranspose(qc::StateCoeffs) = QuCoeffs(conj(qc.arr), flip(qc.tran), flip(qc.conj))
