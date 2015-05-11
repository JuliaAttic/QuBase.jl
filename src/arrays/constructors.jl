############################
# Convenience Constructors #
############################
    statevec(i::Int, fb::FiniteBasis) = QuArray(single_coeff(i, length(fb)), fb)
    statevec(i::Int, lens::Int...=i) = statevec(i, FiniteBasis(lens))

    ###########################
    # Array-like Constructors #
    ###########################
    function Base.zeros(qa::AbstractQuArray)
        fc = zeros(coeffs(qa))
        QAT = similar_type(typeof(qa))
        return QAT(fc, bases(qa))
    end
    function Base.eye(qa::AbstractQuArray)
        fc = eye(coeffs(qa))
        QAT = similar_type(typeof(qa))
        return QAT(fc, bases(qa))
    end

    ####################
    # Helper Functions #
    ####################
    one_at_ind!(arr, i) = setindex!(arr, one(eltype(arr)), i)
    single_coeff(i, lens...) = one_at_ind!(zeros(lens), i)

# Reference :
# Section - 1.1, http://cds.cern.ch/record/331607/files/9708012.pdf
function coherentstatevec_inf(n::Int,alpha::Number)
    s = zeros(typeof(float(alpha)),n)
    s[1] = one(alpha)
    for i in 2:n
        s[i] = alpha/sqrt(i-1)*s[i-1]
    end
    z = QuArray(s)
    return scale!(exp(-abs2(alpha)/2),z)
end

# Reference :
# 2nd Defintion : http://en.wikipedia.org/wiki/Coherent_states#Mathematical_features_of_the_canonical_coherent_states
coherentstatevec(n::Int, alpha::Number) = displaceop(n,alpha)*statevec(1,FiniteBasis(n))

export statevec,
    coherentstatevec
