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
        return spin(@compat(Int(j))) # j = 0.0, 1.0, 2.0...
    elseif isinteger(2*j)
        return HalfSpin(@compat(Int(2*j))) # j = .5, 1.5, 2.5...
    else
        error("Spin must be a multiple of 1/2.")
    end
end

spin_value(s::HalfSpin) = s.val/2

#####################
# Auxiliary Methods #
#####################
spin_coeffs(j, scalar) = [scalar*sqrt(j*(j+1)-(x+1)*x) for x in j-1:-1:-j]

function spinjyx_mat(s::HalfSpin, a, b)
    N = s.val + 1
    j = spin_value(s)
    up_coeffs = spin_coeffs(j, a)
    down_coeffs = b == 1 ? up_coeffs : b * up_coeffs
    return spdiagm((up_coeffs, down_coeffs), (1, -1), N, N)
end

spinjx_mat(j) = spinjyx_mat(spin(j), .5, 1)
spinjy_mat(j) = spinjyx_mat(spin(j), -.5im, -1)

function spinjpm_mat(s::HalfSpin, k)
    N = s.val + 1
    return spdiagm(spin_coeffs(spin_value(s), 1), k, N, N)
end

spinjp_mat(j) = spinjpm_mat(spin(j), 1)
spinjm_mat(j) = spinjpm_mat(spin(j), -1)

function spinjz_mat(s::HalfSpin)
    j = spin_value(s)
    return spdiagm(j:-1:-j)
end

spinjz_mat(j) = spinjz_mat(spin(j))

################################
# Jx, Jy, Jz, Jp, Jm operators #
################################
# These Wikipedia pages give a nice overview of spin/ladder operators:
# http://en.wikipedia.org/wiki/Spin_%28physics%29#Operator
# http://en.wikipedia.org/wiki/Ladder_operator#Angular_momentum
spinjx(j) = QuArray(spinjx_mat(j))
spinjy(j) = QuArray(spinjy_mat(j))
spinjz(j) = QuArray(spinjz_mat(j))
spinjp(j) = QuArray(spinjp_mat(j))
spinjm(j) = QuArray(spinjm_mat(j))

############################
#  Pauli spin 1/2 matrices #
############################
const sigmax = 2.0 * full(spinjx(1/2))
const sigmay = 2.0 * full(spinjy(1/2))
const sigmaz = 2.0 * full(spinjz(1/2))
const sigmap = full(sigmax) + im * full(sigmay)
const sigmam = full(sigmax) - im * full(sigmay)

# TODO : unicode sigmay, sigmaz
const σₓ = sigmax

export spinjx,
    spinjy,
    spinjz,
    spinjp,
    spinjm,
    sigmax,
    sigmay,
    sigmaz,
    sigmap,
    sigmam,
    σₓ
