###################
# AbstractQuArray #
###################
    abstract AbstractQuArray{B<:AbstractBasis,T,N}

    # all subtypes of AbstractQuArray have to implement the following functions
    #   - coefftype
    #   - rawcoeffs
    #   - rawbases
    #   - copy
    #   - similar_type

    coeffs(qarr::AbstractQuArray) = rawcoeffs(qarr)

    # works generally as long as the single index form is defined
    bases(qarr::AbstractQuArray, i) = rawbases(qarr, i)
    bases(qarr::AbstractQuArray) = ntuple(i->bases(qarr, i), ndims(qarr))

    ########################
    # Array-like functions #
    ########################
    Base.eltype{B,T,N}(::AbstractQuArray{B,T,N}) = T
    Base.eltype{B,T,N}(::Type{AbstractQuArray{B,T,N}}) = T

    Base.size(qarr::AbstractQuArray, i...) = size(rawcoeffs(qarr), i...)
    Base.ndims(qarr::AbstractQuArray) = ndims(rawcoeffs(qarr))
    Base.length(qarr::AbstractQuArray) = length(rawcoeffs(qarr))

    Base.getindex(qarr::AbstractQuArray, i...) = getindex(rawcoeffs(qarr), i...)
    Base.setindex!(qarr::AbstractQuArray, i...) = setindex!(rawcoeffs(qarr), i...)

    Base.in(c, qarr::AbstractQuArray) = in(c, rawcoeffs(qarr))

    Base.(:(==))(a::AbstractQuArray, b::AbstractQuArray) = coeffs(a)==coeffs(b) && bases(a)==bases(b)

    typealias AbstractQuVector{B<:AbstractBasis,T} AbstractQuArray{B,T,1}
    typealias AbstractQuMatrix{B<:AbstractBasis,T} AbstractQuArray{B,T,2}

    similar_type{Q<:AbstractQuArray}(::Q) = similar_type(Q)
    similar_type{A<:AbstractQuArray,B<:AbstractQuArray}(::A,::B) = similar_type(promote_type(A,B))

###########
# QuArray #
###########
    type QuArray{B<:AbstractBasis,T,N,A} <: AbstractQuArray{B,T,N}
        coeffs::A
        bases::NTuple{N,B}
        function QuArray(coeffs::AbstractArray{T,N}, bases::NTuple{N,B})
            if checkbases(coeffs, bases)
                new(coeffs, bases)
            else
                error("Coefficient array does not conform to input bases")
            end
        end
    end

    QuArray{B<:AbstractBasis,T,N}(coeffs::AbstractArray{T,N}, bases::NTuple{N,B}) = QuArray{B,T,N,typeof(coeffs)}(coeffs, bases)
    QuArray(coeffs, bases::AbstractBasis...) = QuArray(coeffs, bases)
    QuArray(coeffs) = QuArray(coeffs, map(FiniteBasis, size(coeffs)))

    typealias QuVector{B<:AbstractBasis,T,A} QuArray{B,T,1,A}
    typealias QuMatrix{B<:AbstractBasis,T,A} QuArray{B,T,2,A}

    # returns an appropriate outer constructor for the given (sub)type;
    # the constructor should construct an instance from a coefficient container
    # and a tuple of bases (see  QuArray(coeffs, bases) above)
    similar_type{Q<:QuArray}(::Type{Q}) = QuArray

    ######################
    # Accessor functions #
    ######################
    coefftype{B,T,N,A}(::QuArray{B,T,N,A}) = A
    coefftype{B,T,N,A}(::Type{QuArray{B,T,N,A}}) = A

    rawcoeffs(qarr::QuArray) = qarr.coeffs

    rawbases(qarr::QuArray) = qarr.bases
    rawbases(qarr::QuArray, i) = qarr.bases[i]

    ########################
    # Array-like functions #
    ########################
    Base.copy(qa::QuArray) = QuArray(copy(qa.coeffs), copy(qa.bases))
    Base.eltype{B,T,N,A}(::Type{QuArray{B,T,N,A}}) = T

##############
# CTranspose #
##############
    immutable CTranspose{B<:AbstractBasis,T,N,Q} <: AbstractQuArray{B,T,N}
        qarr::Q
        CTranspose(qarr::AbstractQuVector{B,T}) = new(qarr)
        CTranspose(qarr::AbstractQuMatrix{B,T}) = new(qarr)
        CTranspose(qarr::AbstractQuArray) = error("Conjugate transpose is unsupported for AbstractQuArrays of dimension $N")
    end

    CTranspose{B<:AbstractBasis,T,N}(qa::AbstractQuArray{B,T,N}) = CTranspose{B,T,N,typeof(qa)}(qa)
    CTranspose{B<:AbstractBasis,T,N}(coeffs::AbstractArray{T,N}, bases::NTuple{N,B}) = CTranspose(QuArray(coeffs, bases))

    typealias DualVector{B<:AbstractBasis,T,Q} CTranspose{B,T,1,Q}
    typealias DualMatrix{B<:AbstractBasis,T,Q} CTranspose{B,T,2,Q}

    ######################
    # Accessor functions #
    ######################
    coefftype(ct::CTranspose) = coefftype(ct.qarr)
    coefftype{B,T,N,Q}(::Type{CTranspose{B,T,N,Q}}) = coefftype(Q)

    rawcoeffs(ct::CTranspose) = rawcoeffs(ct.qarr)
    coeffs(ct::CTranspose) = rawcoeffs(ct)'

    revind(len, i) = len - (i-1)

    rawbases(ct::CTranspose, i) = rawbases(ct.qarr, i)
    rawbases(ct::CTranspose) =  rawbases(ct.qarr)

    bases(ct::CTranspose, i) = rawbases(ct, revind(ndims(ct), i))

    qarr_type{B,T,N,Q}(::Type{CTranspose{B,T,N,Q}}) = Q
    qarr_type(ct::CTranspose) = qarr_type(typeof(ct))

    similar_type{QC<:CTranspose}(::Type{QC}) = CTranspose

    ########################
    # Array-like functions #
    ########################
    Base.eltype{B,T,N,Q}(::Type{CTranspose{B,T,N,Q}}) = T

    Base.copy(ct::CTranspose) = CTranspose(copy(ct.qarr))

    Base.ndims(ct::CTranspose) = ndims(ct.qarr)
    Base.length(ct::CTranspose) = length(ct.qarr)

    Base.size(ct::CTranspose) = reverse(size(ct.qarr))
    Base.size(ct::CTranspose, i) = size(ct.qarr, revind(ndims(ct), i))

    Base.getindex(ct::CTranspose, i, j) = getindex(ct.qarr, j, i)'
    Base.getindex(dv::DualVector, i) = getindex(dv.qarr, i)'

    Base.in(c, ct::CTranspose) = in(c', ct.qarr)

    Base.setindex!(ct::CTranspose, x, i, j) = setindex!(ct.qarr, x', j, i)
    Base.setindex!(dv::DualVector, x, i) = setindex!(dv.qarr, x', i)

    Base.ctranspose(qarr::QuArray) = CTranspose(qarr)
    Base.ctranspose(ct::CTranspose) = ct.qarr

    eager_ctranspose(qarr::AbstractQuArray) = similar_type(qarr)(coeffs(qarr)', reverse(bases(qarr)))

################
# LabelQuArray #
################
    typealias LabeledQuArray{B<:LabelBasis,T,N} AbstractQuArray{B,T,N}
    typealias TupleArray{T<:Tuple,N} Array{T,N}

    # redudant definitions added to resolve ambiguity warnings
    Base.getindex{B<:LabelBasis}(larr::DualVector{B}, tups::@compat Union{Tuple,TupleArray}) = getindex(larr, bases(larr,1)[tups])
    Base.getindex{B<:LabelBasis}(larr::CTranspose{B}, tups::@compat(Union{Tuple,TupleArray})...) = getindex(larr, map(getindex, bases(larr), tups)...)
    Base.getindex(larr::LabeledQuArray, tups::@compat(Union{Tuple,TupleArray})...) = getindex(larr, map(getindex, bases(larr), tups)...)

    # redudant definitions added to resolve ambiguity warnings
    Base.setindex!{B<:LabelBasis}(larr::DualVector{B}, x, tups::@compat Union{Tuple,TupleArray}) = setindex(larr, x, bases(larr,1)[tups])
    Base.setindex!{B<:LabelBasis}(larr::CTranspose{B}, x, tups::@compat(Union{Tuple,TupleArray})...) = setindex!(larr, x, map(getindex, bases(larr), tups)...)
    Base.setindex!(larr::LabeledQuArray, x, tups::@compat(Union{Tuple,TupleArray})...) = setindex!(larr, x, map(getindex, bases(larr), tups)...)

########################
# Conversion/Promotion #
########################
    Base.promote_rule{B1,B2,T1,T2,N,A1,A2}(::Type{QuArray{B1,T1,N,A1}}, ::Type{QuArray{B2,T2,N,A2}}) = QuArray{promote_type(B1,B2),promote_type(T1,T2),N,promote_type(A1,A2)}
    Base.promote_rule{C<:CTranspose,B,T,N,A}(::Type{C}, ::Type{QuArray{B,T,N,A}}) = promote_type(qarr_type(C), QuArray{B,T,N,A})
    Base.promote_rule{DV<:DualVector, QV<:QuVector}(::Type{DV}, ::Type{QV}) = typejoin(DV,QV)

    Base.convert{B,T,N,A}(::Type{QuArray{B,T,N,A}}, qa::QuArray{B}) = QuArray(convert(A, coeffs(qa)), bases(qa))
    Base.convert{B,T,N,A}(::Type{QuArray{B,T,N,A}}, qa::QuArray) = QuArray(convert(A, coeffs(qa)), map(i -> convert(B,i), bases(qa)))

    Base.convert{C<:CTranspose,B,T,N,Q<:AbstractQuArray}(::Type{C}, ct::CTranspose{B,T,N,Q}) = CTranspose(convert(qarr_type(C),ct.qarr))
    Base.convert{B,T,N,Q<:AbstractQuArray}(::Type{Q}, ct::CTranspose{B,T,N,Q}) = eager_ctranspose(ct.qarr)
    Base.convert{DV<:DualVector,B,T,QV<:AbstractQuVector}(::Type{DV}, ct::DualVector{B,T,QV}) = CTranspose(convert(qarr_type(DV),ct.qarr))

######################
# Printing Functions #
######################
    Base.summary{B}(qarr::AbstractQuArray{B}) = "$(sizenotation(size(qarr))) $(typerepr(qarr)) in $B"
    Base.summary(larr::LabeledQuArray) = "$(sizenotation(size(larr))) $(typerepr(larr)) in labeled bases $(bases(larr))"

    function Base.show(io::IO, qarr::AbstractQuArray)
        println(io, summary(qarr)*":")
        println(io, "...coefficients: $(coefftype(qarr))")
        print(io, repr(coeffs(qarr)))
    end

####################
# Helper Functions #
####################
    typerepr(qa::AbstractQuArray) = "$(typeof(qa).name)"
    typerepr(::QuVector) = "QuVector"
    typerepr(::QuMatrix) = "QuMatrix"
    typerepr(::DualVector) = "DualVector"
    typerepr(::DualMatrix) = "DualMatrix"
    typerepr(::QuArray) = "QuArray"

    sizenotation(@compat(tup::Tuple{Int})) = "$(first(tup))-element"
    sizenotation(@compat(tup::Tuple{Vararg{Int}})) = reduce(*, map(s->"$(s)x", tup))[1:end-1]

    function checkbases{N}(coeffs, bases::NTuple{N,AbstractBasis})
        if ndims(coeffs) == length(bases)
            return reduce(&, [checkcoeffs(coeffs, i, bases[i]) for i=1:N])
        end
        return false
    end

export QuArray, QuVector, QuMatrix, CTranspose, DualVector, DualMatrix,
    rawcoeffs,
    coeffs,
    coefftype,
    rawbases,
    bases
