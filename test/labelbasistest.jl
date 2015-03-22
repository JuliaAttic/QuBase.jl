m = [1+1im 2+2im 3+3im 4+4im; 
     5+5im 6+6im 7+7im 8+8im; 
     9+9im 10+10im 11+11im 12+12im;
     13+13im 14+14im 15+15im 16+16im]

qm = QuArray(m, LabelBasis(4,1), LabelBasis(2,2))
qm_qm = tensor(qm,qm)

@assert qm[(2,0),(1,0)] == qm[3,2]
@assert qm'[(1,0),(2,0)] == qm[3,2]'
@assert QuBase.bases(qm_qm) == (LabelBasis(4,1,4,1), LabelBasis(2,2,2,2))

@assert qm_qm[(2,0,1,0),(1,1,0,1)] == qm_qm[7,12]
qm_qm[(2,0,1,0),(1,1,0,1)] = 123
@assert qm_qm[(2,0,1,0),(1,1,0,1)] == 123

x2subspace = [filter((x...)->sum(x...)==2, QuBase.bases(qm_qm,2))...]
@assert x2subspace == [(1,1,0,0),(1,0,1,0),(0,1,1,0),(1,0,0,1),(0,1,0,1),(0,0,1,1)]

@assert qm_qm[(3,0,3,0), x2subspace] == [0+416im 0+392im 0+420im 0+420im 0+450im 0+416im]
qm_qm[(3,0,3,0), x2subspace] = [1 2 3 4 5 6]
@assert qm_qm[(3,0,3,0), x2subspace] == [1 2 3 4 5 6]