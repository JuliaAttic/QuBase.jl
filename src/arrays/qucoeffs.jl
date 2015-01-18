import Base: transpose,
    ctranspose,
    size,
    ndims,
    length,
    getindex

############
# QuCoeffs #                 
############
    abstract AbstractQuCoeffs{T,N,A}
    abstract ConjBool{Conj}
    abstract TranBool{Tran}

    type QuCoeffs{Conj,Tran,T,N,A} <: AbstractQuCoeffs{T,N,A}
        arr::A  
        conj::Type{ConjBool{Conj}}
        tran::Type{TranBool{Tran}}
    end

    function QuCoeffs{T,N,Conj,Tran}(arr::AbstractArray{T,N}, 
                                     conj::Type{ConjBool{Conj}}, 
                                     tran::Type{TranBool{Tran}})
        return QuCoeffs{Conj,Tran,T,N,typeof(arr)}(arr, conj, tran)
    end

    QuCoeffs(arr) = QuCoeffs(arr, ConjBool{false}, TranBool{false})
    
    typealias CoeffsVector{Conj,Tran,T,A} QuCoeffs{Conj,Tran,T,1,A}
    typealias CoeffsMatrix{Conj,Tran,T,A} QuCoeffs{Conj,Tran,T,2,A}

    typealias ConjVector{Tran,T,A} CoeffsVector{true,Tran,T,A}
    typealias ConjMatrix{Tran,T,A} CoeffsMatrix{true,Tran,T,A}

    typealias TranVector{Conj,T,A} CoeffsVector{Conj,true,T,A}
    typealias TranMatrix{Conj,T,A} CoeffsMatrix{Conj,true,T,A}

    typealias AdjVector{T,A} CoeffsVector{true,true,T,A}
    typealias AdjMatrix{T,A} CoeffsMatrix{true,true,T,A}

    ########################
    # Array-like Functions #
    ########################
    size(qc::QuCoeffs) = size(qc.arr)
    size(qc::QuCoeffs, i...) = size(qc.arr, i...)
    size(tm::TranMatrix) = reverse(size(tm.arr))
    size(tm::TranMatrix, i) = size(tm)[i]
    ndims(qc::QuCoeffs) = length(size(qc))
    length(qc::QuCoeffs) = prod(size(qc))

    apply_conj(i::Complex, ::Type{ConjBool{true}}) = conj(i)
    apply_conj(i, conj) = i

    getindex(cv::CoeffsVector, i) = apply_conj(cv.arr[i], cv.conj)
    getindex(cv::CoeffsMatrix, i, j) = apply_conj(cv.arr[i,j], cv.conj)
    getindex(cv::TranMatrix, i, j) = apply_conj(cv.arr[j,i], cv.conj)

    #######################
    # Conjugate/Transpose #                 
    #######################
    conj{Conj,Tran}(qc::QuCoeffs{Conj,Tran}) = QuCoeffs(qc.arr, ConjBool{!(Conj)}, TranBool{Tran})
    transpose{Conj,Tran}(qc::QuCoeffs{Conj,Tran}) = QuCoeffs(qc.arr, ConjBool{Conj}, TranBool{!(Tran)})
    ctranspose{Conj,Tran}(qc::QuCoeffs{Conj,Tran}) = QuCoeffs(qc.arr, ConjBool{!(Conj)}, TranBool{!(Tran)})

