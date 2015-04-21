#######################
# Spin Operators Test #
#######################

@assert commute(sigma_x, sigma_y) == 2*im*sigma_z
@assert commute(sigma_y, sigma_z) == 2*im*sigma_x
@assert commute(sigma_z, sigma_x) == 2*im*sigma_y

@assert anticommute(sigma_x, sigma_x) == 2 * sigma_x^2

@assert coeffs(commute(sigma_x, sigma_x)) == spzeros(2,2)

####################################################
# Position, Displacement & Momentum Operators Test #
####################################################

p = positionop(2)
m = momentumop(2)
@assert coeffs(commute(sigma_x, p)) == spzeros(2,2)
@assert coeffs(commute(sigma_y, m)) == spzeros(2,2)

coherentstate_inf = QuBase.coherentstatevec_inf(20,1)
coherentstate = displaceop(20,1)*statevec(1,FiniteBasis(20))
@test_approx_eq_eps coeffs(coherentstate_inf)[1] coeffs(coherentstate)[1] 1e-8

############################
# Squeezing Operators Test #
############################

@assert squeezingop(2,1.0)== QuArray(eye(2))
@test_approx_eq coeffs(squeezingop(lowerop(2), QuArray(eye(2)), 2.0)') coeffs(displaceop(2,1.0))
