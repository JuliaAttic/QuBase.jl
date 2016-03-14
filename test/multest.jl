m = [1+1im 2+2im 3+3im;
     4+4im 5+5im 6+6im;
     7+7im 8+8im 9+9im]

v = [1.0,2.0,3.0]

qm = QuArray(m)
qv = QuArray(v)

@assert coeffs(qm*qm) == m*m
@assert coeffs(qm'*qm) == m'*m
@assert coeffs(qm*qm') == m*m'
@assert coeffs(qm'*qm') == m'*m'
@assert rawcoeffs(qm'*qm') == m*m

@assert coeffs(qm*qv) == m*v
@assert rawcoeffs(qv'*qm') == m*v
@assert coeffs(qv'*qm') == v'*m'
@assert coeffs(qm'*qv) == m'*v
@assert coeffs(qv'*qm) == v'*m

@assert qv'*qv == dot(v,v)
@assert coeffs(qv*qv') == v*v'

@assert qv'*qm*qv == first(v'*m*v)

@assert coeffs(2*im*qv) == 2*im*v
@assert coeffs(qv/2) == v/2
@assert coeffs(scale(2,qv)) == 2*v
@assert coeffs(qm' * -im) == -im * m'
@assert norm(qv) == norm(v)
@test_approx_eq coeffs(normalize(qv)) coeffs(qv)/norm(v)

# a simple test of the `==` operator
@assert qm*3im == scale!(3im, copy(qm))

# tests for Matrix Division
v1 =  [0.5+0*im, 0.+0.*im]
qv1 = normalize!(QuArray(v1))
@assert coeffs(\(sigmax,qv1)) == [0.+0.*im, 1.+0.*im]
@assert coeffs(\(sigmaz, sigmax)) == [0. 1.;-1. 0.]

# Trace
# Pauli matrices sigmax, sigmay, sigmaz are traceless
# Ref : http://en.wikipedia.org/wiki/Pauli_matrices#Algebraic_properties
@assert trace(sigmax) == zero(eltype(sigmax))
@assert trace(sigmax) == trace(sigmaz)

# Dot product
# Being defined both for vectors and dual vectors
# Related Ref : https://github.com/JuliaLang/julia/issues/11064
@assert dot(qv,qv) == dot(qv',qv')
# dot(A,B)=trace(A'*B) for the case of sigmax, sigmay results in no trace.
@assert dot(sigmax, sigmay) == trace(sigmaz)

# Vectorize
@assert vec(qv) == vec(qv')

# Expectation value
@assert expectationvalue(qv1, sigmax) == qv1'*sigmax*qv1
@assert expectationvalue(qm, lowerop(3)) == trace(qm*lowerop(3))
