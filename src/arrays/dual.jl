#####################################
# Conjugate/Transpose Wrapper types #
#####################################
    abstract DualWrapper{B,T,N,A} <: AbstractQuArray{B,T,N}

    immutable Transpose{B,T,N,A} <: DualWrapper{B,T,N,A}
        data::QuArray{B,T,N,A}
        Transpose(data::QuVector{B,T,A}) = new(data)
        Transpose(data::QuMatrix{B,T,A}) = new(data)
        Transpose(data::QuArray) = error("Transposition is unsupported for QuArrays of dimension $N")
    end

    immutable Conjugate{B,T,N,A} <: DualWrapper{B,T,N,A}
        data::QuArray{B,T,N,A}
    end

    immutable CTranspose{B,T,N,A} <: DualWrapper{B,T,N,A}
        data::QuArray{B,T,N,A}
        CTranspose(data::QuVector{B,T,A}) = new(data)
        CTranspose(data::QuMatrix{B,T,A}) = new(data)
        CTranspose(data::QuArray) = error("Conjugate transpose is unsupported for QuArrays of dimension $N")
    end

    Transpose{B,T,N,A}(qa::QuArray{B,T,N,A}) = Transpose{B,T,N,A}(qa)
    Conjugate{B,T,N,A}(qa::QuArray{B,T,N,A}) = Conjugate{B,T,N,A}(qa)
    CTranspose{B,T,N,A}(qa::QuArray{B,T,N,A}) = CTranspose{B,T,N,A}(qa)

    typealias TranVector{B,T,A} Transpose{B,T,1,A}
    typealias TranMatrix{B,T,A} Transpose{B,T,2,A}

    typealias ConjVector{B,T,A} Conjugate{B,T,1,A}
    typealias ConjMatrix{B,T,A} Conjugate{B,T,2,A}

    typealias CTranVector{B,T,A} CTranspose{B,T,1,A}
    typealias CTranMatrix{B,T,A} CTranspose{B,T,2,A}

    ######################
    # Accessor functions #
    ######################
    data(tarr::Transpose) = tarr.data
    data(carr::Conjugate) = carr.data
    data(ctarr::CTranspose) = ctarr.data

    coeffs(dw::DualWrapper) = coeffs(data(dw))
    coefftype{B,T,N,A}(::DualWrapper{B,T,N,A}) = A

    revind(len, i) = len - (i-1)
    bases(tmat::Union(CTranMatrix, TranMatrix), i) = bases(data(tmat), revind(ndims(tmat), i))

    isconj(::Union(Conjugate, CTranspose)) = true
    isconj(x) = false

    istran(::Union(Transpose, CTranspose)) = true
    istran(x) = false

    ########################
    # Array-like functions #
    ########################
    Base.size(dw::DualWrapper, i...) = size(data(dw), i...)
    Base.size(tmat::Union(TranMatrix,CTranMatrix)) = reverse(size(data(tmat)))
    Base.size(tmat::Union(TranMatrix,CTranMatrix), i) = size(tmat, revind(ndims(tmat), i))

    Base.ndims(dw::DualWrapper) = ndims(data(dw))
    Base.length(dw::DualWrapper) = length(data(dw))

    Base.getindex(tarr::Transpose, i, j) = getindex(data(tarr), j, i).'
    Base.getindex(tvec::TranVector, i) = getindex(data(tvec), i).'
    Base.getindex(ctarr::CTranspose, i, j) = getindex(data(ctarr), j, i)'
    Base.getindex(ctvec::CTranVector, i) = getindex(data(ctvec), i)'
    Base.getindex(carr::Conjugate, i...) = conj(getindex(data(carr), i...))

    Base.setindex!(tarr::Transpose, x, i, j) = setindex!(data(tarr), x.', j, i)
    Base.setindex!(tvec::TranVector, x, i) = setindex!(data(tvec), x.', i)
    Base.setindex!(ctarr::CTranspose, x, i, j) = setindex!(data(ctarr), x', j, i)
    Base.setindex!(ctvec::CTranVector, x, i) = setindex!(data(ctvec), x', i)
    Base.setindex!(carr::Conjugate, x, i...) = setindex!(data(carr), conj(x), i...)

    Base.conj(qarr::QuArray) = Conjugate(qarr)
    Base.conj(ctarr::CTranspose) = Transpose(data(ctarr))
    Base.conj(tarr::Transpose) = CTranspose(data(tarr))
    Base.conj(carr::Conjugate) = data(carr)

    Base.transpose(qarr::QuArray) = Transpose(qarr)
    Base.transpose(ctarr::CTranspose) = Conjugate(data(ctarr))
    Base.transpose(tarr::Transpose) = data(tarr)
    Base.transpose(carr::Conjugate) = CTranspose(data(carr))

    Base.ctranspose(qarr::AbstractQuArray) = transpose(conj(qarr))
