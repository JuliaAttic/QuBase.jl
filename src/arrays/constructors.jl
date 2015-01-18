import Base: zeros,
    eye

############################
# Convenience Constructors #
############################
    statevec(i::Int, fb::FiniteBasis) = QuArray(single_coeff(i, length(fb)), fb)
    statevec(i::Int, lens::Int...=i) = statevec(i, FiniteBasis(lens))

    ###########################
    # Array-like Constructors #
    ###########################
    zeros(qa::AbstractQuArray) = QuArray(zeros(coeffs(qa)), bases(qa))
    eye(qa::AbstractQuArray) = QuArray(eye(coeffs(qa)), bases(qa))

    ####################
    # Helper Functions #
    ####################
    one_at_ind!(arr, i) = setindex!(arr, one(eltype(arr)), i)
    single_coeff(i, lens...) = one_at_ind!(zeros(lens), i)

export statevec