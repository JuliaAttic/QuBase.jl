import Base: *, +, -, /

##################
# Multiplication #
##################

# The below needs to be defined more generically;
# the below only defines multiplication on orthonormal,
# N<3 objects

typealias OrthonormalBasis{S<:Orthonormal} AbstractBasis{S}

# bra * ket -> scalar
function *{B<:OrthonormalBasis}(bra::DualVector{B}, ket::AbstractQuVector{B})
    return dot(rawcoeffs(bra), rawcoeffs(ket))
end

# bra * operator -> bra
function *{B<:OrthonormalBasis}(bra::DualVector{B}, op::AbstractQuMatrix{B})
    return QuArray(Ac_mul_B(rawcoeffs(op), rawcoeffs(bra)), bases(op,2))'
end

function *{B<:OrthonormalBasis}(bra::DualVector{B}, op::DualMatrix{B})
    return QuArray(rawcoeffs(op)*rawcoeffs(bra), bases(op,2))'
end

# operator * ket -> ket
function *{B<:OrthonormalBasis}(op::AbstractQuMatrix{B}, ket::AbstractQuVector{B})
    vc = rawcoeffs(op)*rawcoeffs(ket)
    QAT = similar_type(ket)
    return QAT(vc, bases(op,1))
end

function *{B<:OrthonormalBasis}(op::DualMatrix{B}, ket::AbstractQuVector{B})
    vc = Ac_mul_B(rawcoeffs(op), rawcoeffs(ket))
    QAT = similar_type(ket)
    return QAT(vc, bases(op,1))
end

# ket * bra -> operator
function *{B<:OrthonormalBasis}(ket::QuVector{B}, bra::DualVector{B})
    return QuArray(A_mul_Bc(rawcoeffs(ket), rawcoeffs(bra)),
                   bases(ket,1),
                   bases(bra,1))
end

# operator * operator -> operator
*(dm1::DualMatrix, dm2::DualMatrix) = (dm2.qarr*dm1.qarr)'
*{B<:OrthonormalBasis}(dm1::DualMatrix{B}, dm2::DualMatrix{B}) = (dm2.qarr*dm1.qarr)'

function *{B<:OrthonormalBasis}(dm::DualMatrix{B}, qm::AbstractQuMatrix{B})
    mc = Ac_mul_B(rawcoeffs(dm), rawcoeffs(qm))
    QAT = similar_type(qm)
    return QAT(mc, bases(dm,1), bases(qm,2))
end

function *{B<:OrthonormalBasis}(qm::AbstractQuMatrix{B}, dm::DualMatrix{B})
    mc = A_mul_Bc(rawcoeffs(dm), rawcoeffs(qm))
    QAT = similar_type(qm)
    return QAT(mc, bases(qm,1), bases(dm,2))
end

function *{B<:OrthonormalBasis}(qm1::AbstractQuMatrix{B}, qm2::AbstractQuMatrix{B})
    mc = rawcoeffs(qm1)*rawcoeffs(qm2)
    QAT = similar_type(qm1, qm2)
    return QAT(mc, bases(qm1,1), bases(qm2,2))
end

# addition and subtraction
function +{B<:AbstractBasis}(qarr1::AbstractQuArray{B}, qarr2::AbstractQuArray{B})
    if bases(qarr1) == bases(qarr2)
        sc = coeffs(qarr1) + coeffs(qarr2)
        QAT = similar_type(qarr1, qarr2)
        return QAT(sc, bases(qarr1))
    else
        error("Bases not compatible")
    end
end

function -{B<:AbstractBasis}(qarr1::AbstractQuArray{B}, qarr2::AbstractQuArray{B})
    if bases(qarr1) == bases(qarr2)
        sc = coeffs(qarr1) - coeffs(qarr2)
        QAT = similar_type(qarr1, qarr2)
        return QAT(sc, bases(qarr1))
    else
        error("Bases not compatible")
    end
end

+(ct1::CTranspose, ct2::CTranspose) = (ct1.qarr+ct2.qarr)'
-(ct1::CTranspose, ct2::CTranspose) = (ct1.qarr-ct2.qarr)'

# scaling
Base.scale!(num::Number, qarr::AbstractQuArray) = (scale!(num, rawcoeffs(qarr)); return qarr)
Base.scale!(num::Number, ct::CTranspose) = CTranspose(scale!(num', ct.qarr))
Base.scale!(qarr::AbstractQuArray, num::Number) = scale!(num, qarr)

function Base.scale(num::Number, qarr::AbstractQuArray)
    fc = scale(num, rawcoeffs(qarr))
    QAT = similar_type(qarr)
    return QAT(fc, bases(qarr))
end
Base.scale(num::Number, ct::CTranspose) = CTranspose(scale(num', ct.qarr))
Base.scale(qarr::AbstractQuArray, num::Number) = scale(num, qarr)

*(num::Number, qarr::AbstractQuArray) = scale(num, qarr)
*(qarr::AbstractQuArray, num::Number) = scale(qarr, num)
/(qarr::AbstractQuArray, num::Number) = scale(1/num, qarr)

# matrix operations returning a scalar
# normalization
Base.norm(qarr::AbstractQuArray) = vecnorm(rawcoeffs(qarr))

function normalize!(qarr::AbstractQuArray)
    scale!(1/norm(qarr), rawcoeffs(qarr))
    return qarr
end

normalize(qarr::AbstractQuArray) = normalize!(copy(qarr))

# matrix operations returning an array
# sparse to dense
function Base.full(qarr::AbstractQuMatrix)
    fc = full(rawcoeffs(qarr))
    QAT = similar_type(qarr)
    return QAT(fc, bases(qarr))
end
#Base.full(ct::CTranspose) = full(ct.qarr)'

# exponential of dense matrix
function Base.expm(qarr::AbstractQuMatrix)
    fc = expm(full(rawcoeffs(qarr)))
    QAT = similar_type(qarr)
    return QAT(fc, bases(qarr))
end
#Base.expm(ct::CTranspose) = expm(ct.qarr)'

##################
# Tensor Product #
##################

# General tensor product definitions
function tensor{B1,B2,T1,T2,N}(qarr1::AbstractQuArray{B1,T1,N}, qarr2::AbstractQuArray{B2,T2,N})
    tc = kron(coeffs(qarr1), coeffs(qarr2))
    QAT = similar_type(qarr1, qarr2)
    return QAT(tc, map(tensor, bases(qarr1), bases(qarr2)))
end

# defined to resolve ambiguity warnings
tensor(ct1::CTranspose, ct2::CTranspose) = tensor(ct1.qarr, ct2.qarr)'
tensor(ct1::DualVector, ct2::DualVector) = tensor(ct1.qarr, ct2.qarr)'

function tensor(ket::AbstractQuVector, bra::DualVector)
    tc = kron(coeffs(ket), coeffs(bra))
    QAT = similar_type(ket, bra)
    return QAT(tc, bases(ket,1), bases(bra,1))
end

tensor(bra::DualVector, ket::AbstractQuVector) = tensor(ket, bra)

###############
# Commutators #
###############

# (anti)commute two QuMatrices and return result
commute(a::AbstractQuMatrix, b::AbstractQuMatrix) = (a*b) - (b*a)
anticommute(a::AbstractQuMatrix, b::AbstractQuMatrix) = (a*b) + (b*a)

export normalize,
    normalize!,
    tensor,
    commute,
    anticommute
