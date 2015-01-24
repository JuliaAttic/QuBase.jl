import Base: transpose,
    ctranspose,
    conj,
    size,
    ndims,
    length,
    getindex,
    setindex!

############
# QuCoeffs #                 
############
    abstract ConjBool{Conj}
    abstract TranBool{Tran}

    type QuCoeffs{Conj,Tran,T,N,A}
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

    typealias KetCoeffs{T,A} CoeffsVector{false,false,T,A}
    typealias BraCoeffs{T,A} CoeffsVector{true,true,T,A}

    typealias TranCoeffs{Conj,T,N,A} QuCoeffs{Conj,true,T,N,A}
    typealias TranVector{Conj,T,A} TranCoeffs{Conj,T,1,A}
    typealias TranMatrix{Conj,T,A} TranCoeffs{Conj,T,2,A}

    typealias ConjCoeffs{Tran,T,N,A} QuCoeffs{true,Tran,T,N,A}
    typealias ConjVector{Tran,T,A} ConjCoeffs{Tran,T,1,A}
    typealias ConjMatrix{Tran,T,A} ConjCoeffs{Tran,T,2,A}
    
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
    apply_conj(i::Complex, qc::QuCoeffs) = apply_conj(i, qc.conj)
    apply_conj(i, conj) = i

    getindex(cv::CoeffsVector, i) = apply_conj(cv.arr[i], cv)
    getindex(cv::CoeffsMatrix, i, j) = apply_conj(cv.arr[i,j], cv)
    getindex(cv::TranMatrix, i, j) = apply_conj(cv.arr[j,i], cv)

    setindex!(cv::CoeffsVector, x, y) = setindex!(cv, apply_conj(x, cv), y)
    setindex!(cv::CoeffsMatrix, x, y, z) = setindex!(cv, apply_conj(x, cv), y, z)
    setindex!(cv::TranMatrix,  x, y, z) = setindex!(cv, apply_conj(x, cv), z, y)

    #######################
    # Conjugate/Transpose #                 
    #######################
    conj{Conj,Tran}(qc::QuCoeffs{Conj,Tran}) = QuCoeffs(qc.arr, ConjBool{!(Conj)}, TranBool{Tran})
    transpose{Conj,Tran}(qc::QuCoeffs{Conj,Tran}) = QuCoeffs(qc.arr, ConjBool{Conj}, TranBool{!(Tran)})
    ctranspose{Conj,Tran}(qc::QuCoeffs{Conj,Tran}) = QuCoeffs(qc.arr, ConjBool{!(Conj)}, TranBool{!(Tran)})
