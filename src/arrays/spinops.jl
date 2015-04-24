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

function spin_jxy_mat(s::HalfSpin, x_or_y)
    N = s.val + 1
    j = spin_value(s)

    if x_or_y == :x
        scalar = .5
    elseif x_or_y == :y
        scalar = -.5im
    else
        error("spin_jxy_mat can only construct matrices for :x or :y")
    end

    up_coeffs = spin_coeffs(j, scalar)
    down_coeffs = x_or_y == :x ? up_coeffs : -1*up_coeffs
    return spdiagm((up_coeffs, down_coeffs), (1, -1), N, N)
end

function spin_jpm_mat(s::HalfSpin, p_or_m)
    N = s.val + 1
    
    if p_or_m == :p
        k = 1
    elseif p_or_m == :m
        k = -1 
    else
        error("spin_jpm_mat can only construct matrices for :p or :m")
    end
    
    return spdiagm(spin_coeffs(spin_value(s), 1), k, N, N)
end

function spin_jz_mat(s::HalfSpin)
    j = spin_value(s)
    return spdiagm(j:-1:-j)
end

################################
# Jx, Jy, Jz, Jp, Jm operators #
################################
spin_Jx(j) = QuArray(spin_jxy_mat(spin(j), :x))
spin_Jy(j) = QuArray(spin_jxy_mat(spin(j), :y))
spin_Jz(j) = QuArray(spin_jz_mat(spin(j)))
spin_Jp(j) = QuArray(spin_jpm_mat(spin(j), :p))
spin_Jm(j) = QuArray(spin_jpm_mat(spin(j), :m))

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
