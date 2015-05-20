# The purpose of this file is to test 
# the basic (i.e. non-mathematical) functions 
# of values with type Q<:AbstractQuArray.

###########
# QuArray #
###########
qa_coeffs = rand(Complex{Float64},4,8)
qa_basis = tuple(FiniteBasis(2,2), FiniteBasis(4,2))
qa = QuArray(qa_coeffs, qa_basis)

@assert coeffs(qa) == qa_coeffs
@assert rawcoeffs(qa) == qa_coeffs
@assert coefftype(qa) == typeof(qa_coeffs)
@assert QuBase.similar_type(qa) == QuArray
@assert length(qa) == length(qa_coeffs)
@assert size(qa) == size(qa_coeffs)
@assert size(qa,1) == size(qa_coeffs,1)
@assert ndims(qa) == ndims(qa_coeffs)

rand_i = rand(1:size(qa,1))
rand_j = rand(1:size(qa,2))

@assert qa[rand_i, rand_j] == qa_coeffs[rand_i, rand_j]

c = rand(Complex{Float64})
qa[rand_i, rand_j] = c

@assert qa[rand_i, rand_j] == qa_coeffs[rand_i, rand_j] == c
@assert in(qa[rand_i, rand_j], qa)
@assert qa == copy(qa)

@assert bases(qa) == qa_basis
@assert bases(qa,1) == qa_basis[1]
@assert rawbases(qa) == qa_basis
@assert rawbases(qa,2) == qa_basis[2]

mock_type = QuArray{FiniteBasis{Orthonormal}, BigFloat, 2, Matrix{BigFloat}}
prom_T = promote_type(eltype(qa), eltype(mock_type))
prom_A = promote_type(coefftype(qa), coefftype(mock_type))
prom_result = QuArray{FiniteBasis{Orthonormal}, prom_T, 2, prom_A}

@assert promote_type(typeof(qa), mock_type) == prom_result

##############
# CTranspose #
##############
qac = qa'

@assert coeffs(qac) == qa_coeffs'
@assert rawcoeffs(qac) == qa_coeffs
@assert coefftype(qac) == coefftype(qa)
@assert QuBase.similar_type(qac) == CTranspose
@assert length(qac) == length(qa)
@assert size(qac) == reverse(size(qa))
@assert size(qac,1) == size(qa,2)
@assert ndims(qac) == ndims(qa)
@assert qac[rand_j, rand_i] == qa[rand_i, rand_j]'

c = rand(Complex{Float64})
qac[rand_j, rand_i] = c

@assert qac[rand_j, rand_i] == qa[rand_i, rand_j]' == c
@assert in(qac[rand_j, rand_i], qac)

@assert qac == copy(qac)

@assert bases(qac) == reverse(bases(qa))
@assert bases(qac,1) == bases(qa,2)
@assert rawbases(qac) == bases(qa)
@assert rawbases(qac,2) == bases(qa,2)

qac_eager = QuBase.eager_ctranspose(qa)

@assert convert(QuBase.qarr_type(qac), qac) == qac_eager
@assert qac == qac_eager
@assert promote_type(typeof(qac), mock_type) == prom_result
