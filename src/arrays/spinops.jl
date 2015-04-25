#########################
# Spin Operator Methods #
#########################
immutable HalfSpin{T<:Integer}
    val::T
    HalfSpin(val) = val >= 0 ? new(val) : error("Spin must be positive")
end

HalfSpin{T}(val::T) = HalfSpin{T}(val)

spin(s::HalfSpin) = s
spin(j::Integer) = HalfSpin(2*j)
function spin(j)
    if isinteger(j)       
        return spin(int(j)) # j = 0.0, 1.0, 2.0...
    elseif isinteger(2*j) 
        return HalfSpin(int(2*j)) # j = .5, 1.5, 2.5...
    else
        error("Spin must be a multiple of 1/2.")
    end
end

spin_value(s::HalfSpin) = s.val/2

#####################
# Auxiliary Methods #
#####################
spin_coeffs(j, scalar) = [scalar*sqrt(j*(j+1)-(x+1)*x) for x in j-1:-1:-j]

function spin_jyx_mat(s::HalfSpin, a, b)
    N = s.val + 1
    j = spin_value(s)
    up_coeffs = spin_coeffs(j, a)
    down_coeffs = b == 1 ? up_coeffs : b * up_coeffs
    return spdiagm((up_coeffs, down_coeffs), (1, -1), N, N)
end

spin_jx_mat(j) = spin_jyx_mat(spin(j), .5, 1)
spin_jy_mat(j) = spin_jyx_mat(spin(j), -.5im, -1)

function spin_jpm_mat(s::HalfSpin, k)
    N = s.val + 1
    return spdiagm(spin_coeffs(spin_value(s), 1), k, N, N)
end

spin_jp_mat(j) = spin_jpm_mat(spin(j), 1)
spin_jm_mat(j) = spin_jpm_mat(spin(j), -1)

function spin_jz_mat(s::HalfSpin)
    j = spin_value(s)
    return spdiagm(j:-1:-j)
end

spin_jz_mat(j) = spin_jz_mat(spin(j))

################################
# Jx, Jy, Jz, Jp, Jm operators #
################################
# This Wikipedia pages give a nice overview of spin/ladder operators:
# http://en.wikipedia.org/wiki/Spin_%28physics%29#Operator
# http://en.wikipedia.org/wiki/Ladder_operator#Angular_momentum
spin_Jx(j) = QuArray(spin_jx_mat(j))
spin_Jy(j) = QuArray(spin_jy_mat(j))
spin_Jz(j) = QuArray(spin_jz_mat(j))
spin_Jp(j) = QuArray(spin_jp_mat(j))
spin_Jm(j) = QuArray(spin_jm_mat(j))

############################
#  Pauli spin 1/2 matrices #
############################
const sigma_x = 2.0 * spin_Jx(1/2)
const sigma_y = 2.0 * spin_Jy(1/2)
const sigma_z = 2.0 * spin_Jz(1/2)

# TODO : unicode sigma_y, sigma_z
const σₓ = sigma_x

export spin_Jx,
    spin_Jy,
    spin_Jz,
    spin_Jp,
    spin_Jm,
    sigma_x,
    sigma_y,
    sigma_z,
    σₓ
