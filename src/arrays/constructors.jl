import Base: zeros,
    eye

############################
# Convenience Constructors #
############################
    statevec(i::Int, fb::FiniteBasis) = QuArray(single_coeff(i, length(fb)), fb)
    statevec(i::Int, lens::Int...=i) = statevec(i, FiniteBasis(lens))

    zeros(qa::AbstractQuArray) = QuArray(zeros(coeffs(qa)), bases(qa))
    eye(qa::AbstractQuArray) = QuArray(eye(coeffs(qa)), bases(qa))

export statevec