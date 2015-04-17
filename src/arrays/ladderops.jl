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

squeeze_construct(a::AbstractQuMatrix,b::AbstractQuMatrix, z::Number) = scale!(0.5,(scale!(z',a*b)-scale!(z,a'*b')))
squeezingop(a::AbstractQuMatrix, b::AbstractQuMatrix, z::Number) = expm(squeeze_construct(a,b,z))
squeezingop(n::Int,z::Number) = squeezingop(lowerop(n),lowerop(n),z)

displaceop(n::Int,alpha::Number) = expm(scale(alpha,lowerop(n)')-scale(alpha',lowerop(n)))

export raiseop,
    lowerop,
    positionop,
    momentumop,
    squeezingop,
    displaceop
