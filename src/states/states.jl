##############
# QuStateVec #
##############

# pure state (aka state vector) type
type QuStateVec{B,T,A} <: AbstractQuVector{B,T}
    state::QuVector{B,T,A}

    QuStateVec(state::QuVector{B,T,A}) = new(state)
end

QuStateVec{B,T}(coeffs::AbstractVector{T}, basis::B) = QuStateVec(QuArray(coeffs, basis))
QuStateVec{B,T,A}(qa::QuVector{B,T,A}) = QuStateVec{B,T,A}(qa)

# interface functions for AbstractQuArray
coefftype{B,T,A}(::QuStateVec{B,T,A}) = A

rawcoeffs(qa::QuStateVec) = rawcoeffs(qa.state)
rawbases(qa::QuStateVec) = rawbases(qa.state)
rawbases(qa::QuStateVec, i) = rawbases(qa.state, i)

similar_type{Q<:QuStateVec}(::Type{Q}) = QuStateVec

Base.copy(qs::QuStateVec) = QuStateVec(copy(qs.state))
Base.promote_rule{B1,B2,T1,T2,A1,A2}(::Type{QuStateVec{B1,T1,A1}}, ::Type{QuVector{B2,T2,A2}}) = QuArray{promote_type(B1,B2),promote_type(T1,T2),1,promote_type(A1,A2)}


# extract populations (occupation probabilities)
populations( qs::QuStateVec ) = abs2( coeffs(qs) )

# purity of the state vector
purity( qs::QuStateVec ) = norm(qs)^2


###################
# QuDensityMatrix #
###################

# density matrix type
type QuDensityMatrix{B,T,A} <: AbstractQuMatrix{B,T}
    state::QuMatrix{B,T,A}

    QuDensityMatrix(state::QuMatrix{B,T,A}) = new(state)
end

QuDensityMatrix{B,T,A}(qa::QuMatrix{B,T,A}) = QuDensityMatrix{B,T,A}(qa)
QuDensityMatrix{B,T}(coeffs::AbstractMatrix{T}, bases::NTuple{2,B}) = QuDensityMatrix(QuArray(coeffs, bases))
QuDensityMatrix{B,T}(coeffs::AbstractMatrix{T}, basis1::B, basis2::B) = QuDensityMatrix(coeffs, (basis1,basis2))

# TO FIX
#QuDensityMatrix(sv::QuStateVec) = QuDensityMatrix(sv*sv')

# interface functions for AbstractQuArray
coefftype{B,T,A}(::QuDensityMatrix{B,T,A}) = A

rawcoeffs(qa::QuDensityMatrix) = rawcoeffs(qa.state)
rawbases(qa::QuDensityMatrix) = rawbases(qa.state)
rawbases(qa::QuDensityMatrix, i) = rawbases(qa.state, i)

similar_type{Q<:QuDensityMatrix}(::Type{Q}) = QuDensityMatrix

Base.copy(qs::QuDensityMatrix) = QuStateVec(copy(qs.state))
Base.promote_rule{B1,B2,T1,T2,A1,A2}(::Type{QuDensityMatrix{B1,T1,A1}}, ::Type{QuMatrix{B2,T2,A2}}) = QuArray{promote_type(B1,B2),promote_type(T1,T2),1,promote_type(A1,A2)}


# extract populations (occupation probabilities)
populations( qs::QuDensityMatrix ) = diag( coeffs(qs) )

# purity of the state
purity( qs::QuDensityMatrix ) = trace( coeffs(qs)*coeffs(qs) )


export QuStateVec, QuDensityMatrix,
    populations, purity
