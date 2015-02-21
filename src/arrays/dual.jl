#####################################
# Conjugate/Transpose Wrapper types #
#####################################
    abstract DualWrapper{B,T,N,A} <: AbstractQuArray{B,T,N}

    immutable Transpose{B,T,N,A,Q<:AbstractQuArray{B,T,N}} <: DualWrapper{B,T,N,A}
        data::Q  
        Transpose(data::AbstractQuVector{B,T}) = new(data)
        Transpose(data::AbstractQuMatrix{B,T}) = new(data)
    end

    immutable Conjugate{B,T,N,A,Q<:AbstractQuArray{B,T,N}} <: DualWrapper{B,T,N,A}
        data::Q
    end

    Transpose{B,T,N}(qa::AbstractQuArray{B,T,N}) = Transpose{B,T,N,coefftype(qa),typeof(qa)}(qa)
    Conjugate{B,T,N}(qa::AbstractQuArray{B,T,N}) = Conjugate{B,T,N,coefftype(qa),typeof(qa)}(qa)

    typealias TranVector{B,T,A} Transpose{B,T,1,A}
    typealias TranMatrix{B,T,A} Transpose{B,T,2,A}

    typealias ConjVector{B,T,A} Conjugate{B,T,1,A}
    typealias ConjMatrix{B,T,A} Conjugate{B,T,2,A}

    typealias ConjTran{B,T,N,A,Q<:Transpose} Conjugate{B,T,N,A,Q}
    typealias TranConj{B,T,N,A,Q<:Conjugate} Transpose{B,T,N,A,Q}

    typealias DualArray{B,T,N,A} Union(ConjTran{B,T,N,A}, TranConj{B,T,N,A})
    typealias DualVector{B,T,A} DualArray{B,T,1,A}
    typealias DualMatrix{B,T,A} DualArray{B,T,2,A}

    ######################
    # Accessor functions #
    ######################
    data(tarr::Transpose) = tarr.data
    data(carr::Conjugate) = carr.data

    coefftype{B,T,N,A}(::DualWrapper{B,T,N,A}) = A
    coeffs(dw::DualWrapper) = coeffs(data(dw))

    revind(len, i) = len - (i-1)
    getbasis(tmat::TranMatrix, i) = getbasis(data(tmat), revind(ndims(tmat), i))

    isconj(::Union(Conjugate, TranConj)) = true
    isconj(x) = false

    istran(::Union(Transpose, ConjTran)) = true
    istran(x) = false

    ########################
    # Array-like functions #
    ########################
    Base.size(dw::DualWrapper, i...) = size(data(dw), i...)
    Base.size(tmat::TranMatrix) = reverse(size(data(tmat)))
    Base.size(tmat::TranMatrix, i) = size(tmat, revind(ndims(tmat), i))

    Base.ndims(dw::DualWrapper) = ndims(data(dw))
    Base.length(dw::DualWrapper) = length(data(dw))

    Base.getindex(carr::Conjugate, i...) = data(carr)[i...]
    Base.getindex(tarr::Transpose, i, j, k...) = data(tarr)[j,i,k...]
    Base.getindex(tvec::TranVector, i) = data(tvec)[i]

    Base.setindex!(carr::Conjugate, x, i...) = (data(carr)[i...] = conj(x))
    Base.setindex!(tarr::Transpose, x, i, j, k...) = (data(tarr)[j,i,k...] = x)
    Base.setindex!(tvec::TranVector, x, i) = (data(tvec)[i] = x)

    Base.conj(arr::AbstractQuArray) = Conjugate(arr)
    Base.conj(tarr::Transpose) = Transpose(conj(data(tarr)))
    Base.conj(carr::Conjugate) = data(carr)

    Base.transpose(arr::AbstractQuArray) = Transpose(arr)
    Base.transpose(carr::Conjugate) = Conjugate(transpose(data(carr)))
    Base.transpose(tarr::Transpose) = data(tarr)

    Base.ctranspose(arr::AbstractQuArray) = conj(transpose(arr))