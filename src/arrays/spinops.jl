#########################
# Spin Operator Methods #
#########################

immutable HalfSpin{S}
    function HalfSpin()
        S > 0 ? new() : error("Spin must be positive")
    end
end

# The `spin` function converts value to parameter i.e., from the value
# domain to the parameter domain making it type unstable.

spin(j::Integer) = j > 0 ? HalfSpin{2*j}() : error("Spin must be positive")

function spin(j::Float64)
    if (j > 0) && isinteger(j)
        spin(int(j))
    elseif (j > 0) && isinteger(2*j)
        HalfSpin{int(2*j)}()
    else
        error("Spin must be positive and a multiple of 1/2.")
    end
end

spin_value{S}(::HalfSpin{S}) = S/2

#####################
# Auxiliary Methods #
#####################

function mat_coeffs_1(j::Float64)
    m = [j-1:-1:-j]
    N = length(m)+1
    m_coeffs = [sqrt(j*(j+1.0) - (x+1.0)*x) for x in m]
    return spdiagm(m_coeffs,1,N,N)
end

function mat_coeffs_2(j::Float64)
    m = [j:-1:-j]
    N = length(m)
    return spdiagm(m,0,N,N)
end

function spin_h1(k::HalfSpin)
    j = spin_value(k)
    return mat_coeffs_1(j)
end

function spin_h2(k::HalfSpin)
    j = spin_value(k)
    return mat_coeffs_2(j)
end

################################
# Jx, Jy, Jz, Jp, Jm operators #
################################

spin_Jx(j) = QuArray(0.5*(spin_h1(spin(j))+spin_h1(spin(j))'))
spin_Jy(j) = QuArray(-0.5*im*(spin_h1(spin(j))-spin_h1(spin(j))'))
spin_Jz(j) = QuArray(spin_h2(spin(j)))
spin_Jp(j) = QuArray(spin_h1(spin(j)))
spin_Jm(j) = QuArray(spin_h1(spin(j))')

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
