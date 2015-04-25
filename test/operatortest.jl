#######################
# Spin Operators Test #
#######################

# the Pauli spin matrices sigma_x, sigma_y and sigma_z obey the commutation
# relations $ [\sigma_a, \sigma_b]  = 2i\epsion_{abc} \sigma_c $, where
# $\epsion_{abc}$ is the Levi-Civita symbol.
# Hence,
@assert commute(sigma_x, sigma_y) == 2*im*sigma_z
@assert commute(sigma_y, sigma_z) == 2*im*sigma_x
@assert commute(sigma_z, sigma_x) == 2*im*sigma_y

@assert coeffs(commute(sigma_x, sigma_x)) == spzeros(2,2)

# Moreover, the sigmas fulfill anticommutation relations
# ${\sigma_a, \sigma_b} = 2\delta_{a,b} I$ wit $\delta_{a,b}$
# being the Kronecker delta and $I$ the 2x2 identity matrix
# Thus,
@assert anticommute(sigma_x, sigma_x) == 2 * QuArray(speye(2))

@assert spin_Jp(3.5) == spin_Jx(3.5) + (im * spin_Jy(3.5))
@assert spin_Jm(3.5) == spin_Jx(3.5) - (im * spin_Jy(3.5))

####################################################
# Position, Displacement & Momentum Operators Test #
####################################################

# position and momentum operators obey the commutation relation
# $[x,p] = i I$, where $I$ is the identity operator
# However, here we just use the fact that for n=2, x~sigma_x and
# p ~ sigma_y
p = positionop(2)
m = momentumop(2)
@assert coeffs(commute(sigma_x, p)) == spzeros(2,2)
@assert coeffs(commute(sigma_y, m)) == spzeros(2,2)

# test coefficients computed by coherentstatevec vs 
# Fock basis coefficients
coherentstate_inf = QuBase.coherentstatevec_inf(20,1.0)
coherentstate = coherentstatevec(20, 1.0)
@test_approx_eq_eps coherentstate_inf[1] coherentstate[1] 1e-8


############################
# Squeezing Operators Test #
############################

@assert squeezingop(2,1.0)== QuArray(eye(2))
@test_approx_eq coeffs(squeezingop(lowerop(2), QuArray(eye(2)), -2.0)) coeffs(displaceop(2,1.0))
