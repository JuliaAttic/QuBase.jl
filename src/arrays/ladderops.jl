###########################
# Ladder Operator Methods #
###########################
    # n specifies which particle (a.k.a tensor product
    # factor) the operator acts on
    raiseop(b::AbstractBasis, n=1) = QuArray(raisematrix(size(b), n), b, b)
    raiseop(lens, n=1) = raiseop(FiniteBasis(lens), n)

    lowerop(b::AbstractBasis, n=1) = QuArray(lowermatrix(size(b), n), b, b)
    lowerop(lens, n=1) = lowerop(FiniteBasis(lens), n)

    ##########################
    # Helper Functions/Types #
    ##########################
    raisematrix(lens, particle) = eye_sandwich(lens, particle, raise_single(lens[particle]))
    lowermatrix(lens, particle) = eye_sandwich(lens, particle, lower_single(lens[particle]))

    laddercoeffs(n) = sqrt(linspace(1, n, n))
    lower_single(n) = sparse([1:n-1], [2:n], laddercoeffs(n-1), n, n)
    raise_single(n) = sparse([2:n], [1:n-1], laddercoeffs(n-1), n, n)

    before_eye(lens, pivot) = speye(prod(lens[[1:pivot-1]]))
    after_eye(lens, pivot) = speye(prod(lens[[pivot+1:length(lens)]]))
    function eye_sandwich(lens, pivot, op)
        return kron(before_eye(lens, pivot),
                    op,
                    after_eye(lens, pivot))
    end

function positionop(n::Int)
    cop = raiseop(n)
    return scale!(1/sqrt(2.),cop+cop')
end

function momentumop(n::Int)
    cop = raiseop(n)
    return scale(im/sqrt(2.), cop-cop')
end

function spin_h1(j)
    m = reverse([-j:j-1])
    N = length(m)+1
    m_coeffs = [sqrt(j*(j+1.0) - (x+1.0)*x) for x in m]
    return spdiagm(m_coeffs,1,N,N)
end

function spin_h2(j)
    m = reverse([-j:j])
    N = length(m)
    return spdiagm(m,0,N,N)
end

function spin_jx(j)
    return QuArray(0.5*(spin_h1(j)+spin_h1(j)'))
end

function spin_jy(j)
    return QuArray(-0.5*im*(spin_h1(j)-spin_h1(j)'))
end

function spin_jz(j)
    return QuArray(spin_h2(j))
end

function spin_jp(j)
    return QuArray(spin_h1(j))
end

function spin_jm(j)
    return QuArray(spin_h1(j)')
end

export raiseop,
    lowerop,
    positionop,
    momentumop,
    spin_jx,
    spin_jy,
    spin_jz,
    spin_jp,
    spin_jm

