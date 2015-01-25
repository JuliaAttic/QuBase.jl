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
    abstract DualBool{D}

    type QuCoeffs{D,T,N,A}
        arr::A  
        dual::Type{DualBool{D}}
        function QuCoeffs(arr::AbstractArray{T,N}, dual::Type{DualBool{D}})
            return new(arr, dual)
        end
    end

    function QuCoeffs{T,N,D}(arr::AbstractArray{T,N}, dual::Type{DualBool{D}})
        return QuCoeffs{D,T,N,typeof(arr)}(arr, dual)
    end

    QuCoeffs(arr::AbstractArray) = QuCoeffs(arr, DualBool{false})

    typealias CoeffsVector{D,T,A} QuCoeffs{D,T,1,A}
    typealias CoeffsMatrix{D,T,A} QuCoeffs{D,T,2,A}

    typealias KetCoeffs{T,A} CoeffsVector{true,T,A}
    typealias BraCoeffs{T,A} CoeffsVector{false,T,A}

    ########################
    # Array-like Functions #
    ########################
    size(cm::CoeffsMatrix{true}) = reverse(size(cm.arr))
    size(cm::CoeffsMatrix{true}, i) = size(cm)[i]
    size(qc::QuCoeffs) = size(qc.arr)
    size(qc::QuCoeffs, i...) = size(qc.arr, i...)

    ndims(qc::QuCoeffs) = length(size(qc))
    length(qc::QuCoeffs) = prod(size(qc))

    apply_conj(i, ::Type{DualBool{true}}) = conj(i)
    apply_conj(i, ::Type{DualBool{false}}) = i

    getindex(cv::CoeffsVector, i) = apply_conj(cv.arr[i], cv.dual)
    getindex(cm::CoeffsMatrix, i::Number, j::Number) = apply_conj(cm.arr[i,j], cm.dual)
    getindex(cm::CoeffsMatrix{true}, i::Number, j::Number) = apply_conj(cm.arr[j,i], cm.dual)

    setindex!(cv::CoeffsVector, x, y) = setindex!(cv.arr, apply_conj(x, cv.dual), y)
    setindex!(cm::CoeffsMatrix, x, y::Number, z::Number) = setindex!(cm.arr, apply_conj(x, cm.dual), y, z)
    setindex!(cm::CoeffsMatrix{true}, x, y::Number, z::Number) = setindex!(cm.arr, apply_conj(x, cm.dual), z, y)

    #######################
    # Conjugate/Transpose #
    #######################
    conj(qc::QuCoeffs) = QuCoeffs(conj(qc.arr), qc.dual)
    transpose(qc::QuCoeffs) = QuCoeffs(transpose(qc.arr), qc.dual)
    ctranspose(qc::QuCoeffs) = QuCoeffs(ctranspose(qc.arr), qc.dual)

    dual{D}(qc::QuCoeffs{D}) = QuCoeffs(qc.arr, DualBool{!D})
