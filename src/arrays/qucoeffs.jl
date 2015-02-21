############
# QuCoeffs #
############
    abstract BoolVal{b}

    # enumerate cases for type-stable flipping
    flip(::Type{BoolVal{false}}) = BoolVal{true}
    flip(::Type{BoolVal{true}}) = BoolVal{false}

    immutable QuCoeffs{Tran,Conj,N,T,A}
        arr::A  
    end

    function QuCoeffs{Conj,Tran,N,T}(arr::AbstractArray{T,N}, 
                      ::Type{BoolVal{Tran}},
                      ::Type{BoolVal{Conj}})
        return QuCoeffs{Tran,Conj,N,T,typeof(arr)}(arr)
    end

    QuCoeffs(arr::AbstractArray) = QuCoeffs(arr, BoolVal{false}, BoolVal{false})

    typealias VecCoeffs{Tran,Conj,T,A} QuCoeffs{Tran,Conj,1,T,A}
    typealias MatCoeffs{Tran,Conj,T,A} QuCoeffs{Tran,Conj,2,T,A}

    typealias ColCoeffs{Conj,T,A} VecCoeffs{false,Conj,T,A}
    typealias RowCoeffs{Conj,T,A} VecCoeffs{true,Conj,T,A}

    typealias KetCoeffs{T,A} ColCoeffs{false,T,A}
    typealias BraCoeffs{T,A} RowCoeffs{true,T,A}

    conjbool{Tran,Conj}(::QuCoeffs{Tran,Conj}) = BoolVal{Conj}
    tranbool{Tran}(::QuCoeffs{Tran}) = BoolVal{Tran}

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
    Base.conj(qc::QuCoeffs) = QuCoeffs(conj(qc.arr), tranbool(qc), flip(conjbool(qc)))
    
    Base.transpose(qc::QuCoeffs) = QuCoeffs(transpose(qc.arr), flip(tranbool(qc)), conjbool(qc))
    Base.transpose(qc::VecCoeffs) = QuCoeffs(copy(qc.arr), flip(tranbool(qc)), conjbool(qc))
 
    Base.ctranspose(qc::QuCoeffs) = QuCoeffs(ctranspose(qc.arr), flip(tranbool(qc)), flip(conjbool(qc)))
    Base.ctranspose(qc::VecCoeffs) = QuCoeffs(conj(qc.arr), flip(tranbool(qc)), flip(conjbool(qc)))
