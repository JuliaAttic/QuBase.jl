using QuBase.bases

m = [1+1im 2+2im 3+3im 4+4im; 
     5+5im 6+6im 7+7im 8+8im; 
     9+9im 10+10im 11+11im 12+12im;
     13+13im 14+14im 15+15im 16+16im]

v = [1+1im, 2+2im, 3+3im, 4+4im]

qm = QuArray(m, FiniteBasis(4,1), FiniteBasis(2,2))
qv = QuArray(v)

# operator ⊗ operator
qm_qm = tensor(qm, qm)
@assert rawcoeffs(qm_qm) == kron(m, m)
@assert bases(qm_qm) == (tensor(bases(qm, 1), bases(qm, 1)), tensor(bases(qm, 2), bases(qm, 2)))

qmc_qm = tensor(qm', qm)
@assert rawcoeffs(qmc_qm) == kron(m', m)
@assert bases(qmc_qm) == (tensor(bases(qm, 2), bases(qm, 1)), tensor(bases(qm, 1), bases(qm, 2)))

qm_qmc = tensor(qm, qm')
@assert rawcoeffs(qm_qmc) == kron(m, m')
@assert bases(qm_qmc) == (tensor(bases(qm, 1), bases(qm, 2)), tensor(bases(qm, 2), bases(qm, 1)))

qmc_qmc = tensor(qm', qm')
@assert rawcoeffs(qmc_qmc) == kron(m, m)
@assert bases(qmc_qmc) == (tensor(bases(qm, 2), bases(qm, 2)), tensor(bases(qm, 1), bases(qm, 1)))

# state ⊗ state
qv_qv = tensor(qv, qv)
@assert rawcoeffs(qv_qv) == kron(v, v)
@assert bases(qv_qv) == (tensor(bases(qv, 1), bases(qv, 1)),)

qv_qvc = tensor(qv, qv')
@assert rawcoeffs(qv_qvc) == kron(v, v')
@assert bases(qv_qvc) == (bases(qv, 1), bases(qv, 1))

qvc_qv = tensor(qv', qv)
@assert rawcoeffs(qvc_qv) == kron(v', v)
@assert bases(qvc_qv) == (bases(qv, 1), bases(qv, 1))

qvc_qvc = tensor(qv', qv')
@assert rawcoeffs(qvc_qvc) == kron(v, v)
@assert bases(qvc_qvc) == (tensor(bases(qv, 1), bases(qv, 1)),)

