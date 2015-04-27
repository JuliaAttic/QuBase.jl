#######################
# Spin Operators Test #
#######################

# the Pauli spin matrices sigmax, sigmay and sigmaz obey the commutation
# relations $ [\sigma_a, \sigma_b]  = 2i\epsion_{abc} \sigma_c $, where
# $\epsion_{abc}$ is the Levi-Civita symbol.
# Hence,
@assert commute(sigmax, sigmay) == 2*im*sigmaz
@assert commute(sigmay, sigmaz) == 2*im*sigmax
@assert commute(sigmaz, sigmax) == 2*im*sigmay

@assert coeffs(commute(sigmax, sigmax)) == spzeros(2,2)

# Moreover, the sigmas fulfill anticommutation relations
# ${\sigma_a, \sigma_b} = 2\delta_{a,b} I$ with $\delta_{a,b}$
# being the Kronecker delta and $I$ the 2x2 identity matrix
# Thus,
@assert anticommute(sigmax, sigmax) == 2 * QuArray(speye(2))

@assert spinjp(3.5) == spinjx(3.5) + (im * spinjy(3.5))
@assert spinjm(3.5) == spinjx(3.5) - (im * spinjy(3.5))

####################################################
# Position, Displacement & Momentum Operators Test #
####################################################

# position and momentum operators obey the commutation relation
# $[x,p] = i I$, where $I$ is the identity operator
# However, here we just use the fact that for n=2, x~sigmax and
# p ~ sigmay
p = positionop(2)
m = momentumop(2)
@assert coeffs(commute(sigmax, p)) == spzeros(2,2)
@assert coeffs(commute(sigmay, m)) == spzeros(2,2)

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
