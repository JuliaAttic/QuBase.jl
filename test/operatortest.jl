#######################
# Spin Operators Test #
#######################

@assert commutator(sigma_x, sigma_y) == 2*im*sigma_z
@assert commutator(sigma_y, sigma_z) == 2*im*sigma_x
@assert commutator(sigma_z, sigma_x) == 2*im*sigma_y
@assert coeffs(commutator(sigma_x, sigma_x)) == spzeros(2,2)

######################################
# Position & Momentum Operators Test #
######################################

p = positionop(2)
m = momentumop(2)
@assert coeffs(commutator(sigma_x, p)) == spzeros(2,2)
@assert coeffs(commutator(sigma_y, m)) == spzeros(2,2)
