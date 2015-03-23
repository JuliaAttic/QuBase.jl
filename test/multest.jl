m = [1+1im 2+2im 3+3im;
     4+4im 5+5im 6+6im;
     7+7im 8+8im 9+9im]

v = [1.0,2.0,3.0]

qm = QuArray(m)
qv = QuArray(v)

@assert rawcoeffs(qm*qm) == m*m
@assert rawcoeffs(qm'*qm') == m*m
@assert rawcoeffs(qm'*qm) == m'*m
@assert rawcoeffs(qm*qm') == m*m'

@assert rawcoeffs(qm*qv) == m*v
@assert rawcoeffs(qv'*qm') == m*v
@assert rawcoeffs(qm'*qv) == m'*v
@assert rawcoeffs(qv'*qm) == m'*v

@assert qv'*qv == dot(v,v)
@assert rawcoeffs(qv*qv') == v*v'

@assert qv'*qm*qv == first(v'*m*v)

@assert rawcoeffs(2*im*qv) == 2*im*v
@assert rawcoeffs(qv/2) == v/2
@assert rawcoeffs(scale(2,qv)) == 2*v
@assert coeffs(qm' * -im) == -im * m'
@assert norm(qv) == norm(v)
@test_approx_eq rawcoeffs(normalize(qv)) rawcoeffs(qv)/norm(v)
