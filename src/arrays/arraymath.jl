import Base: *

##################
# Multiplication #
##################

# The below needs to be defined more generically;
# the below only defines multiplication on orthonormal,
# N<3 objects

typealias OrthonormalBasis{S<:Orthonormal} AbstractBasis{S}

# bra * ket -> scalar
function *{B<:OrthonormalBasis}(bra::DualVector{B}, ket::QuVector{B})
    return dot(rawcoeffs(bra), rawcoeffs(ket))
end

# bra * operator -> bra
function *{B<:OrthonormalBasis}(bra::DualVector{B}, op::QuMatrix{B})
    return QuArray(Ac_mul_B(rawcoeffs(op), rawcoeffs(bra)), bases(op,2))'
end

function *{B<:OrthonormalBasis}(bra::DualVector{B}, op::DualMatrix{B})
    return QuArray(rawcoeffs(op)*rawcoeffs(bra), bases(op,2))'
end

# operator * ket -> ket
function *{B<:OrthonormalBasis}(op::QuMatrix{B}, ket::QuVector{B})
    return QuArray(rawcoeffs(op)*rawcoeffs(ket), bases(op,1))
end

function *{B<:OrthonormalBasis}(op::DualMatrix{B}, ket::QuVector{B})
    return QuArray(Ac_mul_B(rawcoeffs(op), rawcoeffs(ket)), bases(op,1))
end

# ket * bra -> operator
function *{B<:OrthonormalBasis}(ket::QuVector{B}, bra::DualVector{B})
    return QuArray(A_mul_Bc(rawcoeffs(ket), rawcoeffs(bra)),
                   bases(ket,1),
                   bases(bra,1))
end

# operator * operator -> operator
function *{B<:OrthonormalBasis}(qm::QuMatrix{B}, dm::DualMatrix{B})
    return QuArray(A_mul_Bc(rawcoeffs(qm), rawcoeffs(dm)),
                   bases(qm,1),
                   bases(dm,2))
end

function *{B<:OrthonormalBasis}(dm::DualMatrix{B}, qm::QuMatrix{B})
    return QuArray(Ac_mul_B(rawcoeffs(dm), rawcoeffs(qm)),
                   bases(dm,1),
                   bases(qm,2))
end

function *{B<:OrthonormalBasis}(qm1::QuMatrix{B}, qm2::QuMatrix{B})
    return QuArray(rawcoeffs(qm1)*rawcoeffs(qm2),
                   bases(qm1,1),
                   bases(qm2,2))
end

*(dm1::DualMatrix, dm2::DualMatrix) = (dm1.qarr*dm2.qarr)'

# scaling
Base.scale!(num::Number, qarr::QuArray) = (scale!(num, rawcoeffs(qarr)); return qarr)
Base.scale!(num::Number, ct::CTranspose) = CTranspose(scale!(num', ct.qarr))
Base.scale!(qarr::Union(QuArray,CTranspose), num::Number) = scale!(num, qarr)

Base.scale(num::Number, qarr::QuArray) = QuArray(scale(num, rawcoeffs(qarr)), rawbases(qarr))
Base.scale(num::Number, ct::CTranspose) = CTranspose(scale(num', ct.qarr))
Base.scale(qarr::Union(QuArray,CTranspose), num::Number) = scale(num, qarr)

*(num::Number, qarr::Union(QuArray,CTranspose)) = scale(num, qarr)
*(qarr::Union(QuArray,CTranspose), num::Number) = scale(qarr, num)
/(qarr::Union(QuArray,CTranspose), num::Number) = scale(1/num, qarr)

# normalization
Base.norm(qarr::Union(QuArray,CTranspose)) = vecnorm(rawcoeffs(qarr))

function normalize!(qarr::Union(QuArray,CTranspose))
    scale!(1/norm(qarr), rawcoeffs(qarr))
    return qarr
end

normalize(qarr::Union(QuArray,CTranspose)) = normalize!(copy(qarr))

##################
# Tensor Product #
##################

# General tensor product definitions for orthonormal bases
function tensor{B<:OrthonormalBasis,T1,T2,N}(qarr1::AbstractQuArray{B,T1,N}, qarr2::AbstractQuArray{B,T2,N})
    return QuArray(kron(coeffs(qarr1), coeffs(qarr2)), 
                   map(tensor, bases(qarr1), bases(qarr2)))
end

function tensor{B<:OrthonormalBasis,T1,T2,N}(qarr1::CTranspose{B,T1,N}, qarr2::CTranspose{B,T2,N})
    return QuArray(kron(rawcoeffs(qarr1), rawcoeffs(qarr2)), 
                   map(tensor, rawbases(qarr1), rawbases(qarr2)))'
end

function tensor{B<:OrthonormalBasis}(ket::QuVector{B}, bra::DualVector{B})
    return QuArray(kron(coeffs(ket), coeffs(bra)), 
                   bases(ket,1), 
                   bases(bra,1))
end

tensor{B<:OrthonormalBasis}(bra::DualVector{B}, ket::QuVector{B}) = tensor(ket, bra)

export normalize,
    normalize!
