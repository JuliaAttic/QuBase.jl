import Base: *

##################
# Multiplication #
##################
# bra * ket -> scalar
function *{B}(bra::DualVector{B}, ket::QuVector{B})
    return dot(rawcoeffs(bra), rawcoeffs(ket))
end

# bra * operator -> bra
function *{B}(bra::DualVector{B}, op::QuMatrix{B}) 
    return QuArray(Ac_mul_B(rawcoeffs(op), rawcoeffs(bra)), bases(op,2))'
end

function *{B}(bra::DualVector{B}, op::DualMatrix{B}) 
    return QuArray(rawcoeffs(op)*rawcoeffs(bra), bases(op,2))'
end

# operator * ket -> ket
function *{B}(op::QuMatrix{B}, ket::QuVector{B}) 
    return QuArray(rawcoeffs(op)*rawcoeffs(ket), bases(op,1))
end

function *{B}(op::DualMatrix{B}, ket::QuVector{B})
    return QuArray(Ac_mul_B(rawcoeffs(op), rawcoeffs(ket)), bases(op,1))
end

# ket * bra -> operator
function *{B}(ket::QuVector{B}, bra::DualVector{B}) 
    return QuArray(A_mul_Bc(rawcoeffs(ket), rawcoeffs(bra)), 
                   bases(ket,1), 
                   bases(bra,1))
end

# operator * operator -> operator
function *{B}(qm::QuMatrix{B}, dm::DualMatrix{B})
    return QuArray(A_mul_Bc(rawcoeffs(qm), rawcoeffs(dm)), 
                   bases(qm,1), 
                   bases(dm,2))
end

function *{B}(dm::DualMatrix{B}, qm::QuMatrix{B})
    return QuArray(Ac_mul_B(rawcoeffs(dm), rawcoeffs(qm)), 
                   bases(dm,1), 
                   bases(qm,2))
end

function *{B}(qm1::QuMatrix{B}, qm2::QuMatrix{B}) 
    return QuArray(rawcoeffs(qm1)*rawcoeffs(qm2), 
                   bases(qm1,1), 
                   bases(qm2,2))
end

*(dm1::DualMatrix, dm2::DualMatrix) = (dm1.qarr*dm2.qarr)'
