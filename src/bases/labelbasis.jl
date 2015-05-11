##############
# LabelBasis #
##############
    # A LabelBasis uses precomputed values to generate labels 
    # for given indices in the basis, or vice versa (an index 
    # for a given state in the basis).
    #
    # For example:
    #
    #   julia> b = LabelBasis(2,2,2)
    #   LabelBasis{Orthonormal,3}(0:1,0:1,0:1)
    #
    #   julia> collect(b)
    #   8-element Array{Int64,Int64,Int64),1}:
    #    (0,0,0)
    #    (1,0,0)
    #    (0,1,0)
    #    (1,1,0)
    #    (0,0,1)
    #    (1,0,1)
    #    (0,1,1)
    #    (1,1,1)
    #
    #   julia> b[6]
    #   (1,0,1)
    #
    #   julia> b[(1,0,1)]
    #   6
    # 
    # Because the labels are generated rather than actually 
    # stored, one can represent very large bases without 
    # any storage overhead:
    #
    #   julia> b = LabelBasis(221,135,31,42,321,3)
    #   LabelBasis{Orthonormal,6}(0:220,0:134,0:30,0:41,0:320,0:2)
    #   
    #   julia> length(b)
    #   37407898710
    #   
    #   julia> last(b)
    #   (220,134,30,41,320,2)
    #   
    #   julia> b[34234134]
    #   (128,60,0,37,0,0)
    #   
    #   julia> b[(128,60,0,37,0,0)]
    #   34234134
    # 
    # Arbitrary numeric ranges are supported for labels:
    #   
    #   julia> b = LabelBasis(0.0:0.1:0.2, 4:7)
    #   LabelBasis{Orthonormal,2}(0.0:0.1:0.2,4:7)     
    #       
    #   julia> collect(b)
    #   12-element Array{(Float64,Int64),1}:
    #    (0.0,4)
    #    (0.1,4)
    #    (0.2,4)
    #    (0.0,5)
    #    (0.1,5)
    #    (0.2,5)
    #    (0.0,6)
    #    (0.1,6)
    #    (0.2,6)
    #    (0.0,7)
    #    (0.1,7)
    #    (0.2,7)        
    #
    #   julia> b[b[(0.0, 7)]] == (0.0, 7)
    #   true
    #
    # Since LabelBasis is iterable, you can apply filters. The below selects
    # the labels for the X=2 subspace from the given basis:
    #   
    #   julia> collect(filter((x...)->sum(x...)==2, LabelBasis(4,4,4)))
    #   6-element Array{Any,1}:
    #    (2,0,0)
    #    (1,1,0)
    #    (0,2,0)
    #    (1,0,1)
    #    (0,1,1)
    #    (0,0,2)
    
    immutable LabelBasis{S<:AbstractStructure,N} <: AbstractFiniteBasis{S}
        ranges::NTuple{N,Range}
        denoms::NTuple{N,Float64}
        LabelBasis(ranges, denoms, ::Type{BypassFlag}) = new(ranges, denoms)
        function LabelBasis(ranges::NTuple{N,Range})
            # reverse is done to match cartesianmap order
            return LabelBasis{S,N}(ranges, precompute_denoms(reverse(map(length,ranges))), BypassFlag) 
        end

    end

    #resolves ambiguity warnings
    LabelBasis{S<:AbstractStructure}(@compat(::Tuple{}), ::Type{S}=Orthonormal) = error("LabelBasis requires a tuple of ranges as a constructor argument") 
    
    LabelBasis{S<:AbstractStructure,N}(lens::NTuple{N,Range}, ::Type{S}=Orthonormal) = LabelBasis{S,N}(lens)
    LabelBasis{S<:AbstractStructure}(@compat(lens::Tuple{Vararg{Int}}), ::Type{S}=Orthonormal) = LabelBasis(map(n->zero(eltype(n)):(n-1), lens), S)
    LabelBasis(lens...) = LabelBasis(lens)

    Base.convert{A,B,N}(::Type{LabelBasis{A,N}}, b::LabelBasis{B,N}) = LabelBasis{A,N}(b.ranges, b.denoms, BypassFlag)
    Base.convert{A,B,N}(::Type{FiniteBasis{A,N}}, b::LabelBasis{B,N}) = FiniteBasis{A,N}(size(b))

    Base.copy{S,N}(b::LabelBasis{S,N}) = LabelBasis{S,N}(copy(b.ranges), copy(b.denoms), BypassFlag)

    ####################
    # Helper Functions #
    ####################
    # This function precomputes the 
    # denominators for each factor of 
    # the cartesian product.
    #
    # This site offers a thorough 
    # explanation of this method:
    # http://phrogz.net/lazy-cartesian-product
    #
    function precompute_denoms(lens)
        # storing as Floats avoids number precision issues for 
        # outrageously large bases
        total_divisor = prod(float,lens) 
        function get_denom(i)
            total_divisor = div(total_divisor, i)
            return max(1.0, total_divisor)
        end
        # reverse is done to match cartesianmap order
        return reverse(map(get_denom, lens))
    end

    ######################
    # Property Functions #
    ######################
    ranges(b::LabelBasis) = b.ranges
    ranges(b::LabelBasis, i) = b.ranges[i]

    Base.size(b::LabelBasis) = map(length, ranges(b))
    Base.size(b::LabelBasis, i) = length(ranges(b, i))
    Base.length(b::LabelBasis) = prod(length, ranges(b))

    structure{S}(::Type{LabelBasis{S}}) = S
    structure{S,N}(::Type{LabelBasis{S,N}}) = S

    nfactors{S,N}(::LabelBasis{S,N}) = N
    checkcoeffs(coeffs, dim, b::LabelBasis) = size(coeffs, dim) == length(b)

    ######################
    # Accessor Functions #
    ######################
    ind_value(n, range, denom, modulus) = range[@compat(Int(div(n, denom) % modulus))+1]
    tuple_at_ind(b::LabelBasis, i) = ntuple(nfactors(b), x->ind_value(i-1, ranges(b,x), b.denoms[x], size(b,x)))
    pos_in_range(r::Range, i) = i in r ? (i-first(r))/step(r) : throw(BoundsError())
    getpos(b::LabelBasis, label) = @compat Int(sum(map(*, map(pos_in_range, ranges(b), label), b.denoms))+1)

    Base.in(label, b::LabelBasis) = reduce(&, map(in, label, ranges(b)))

    Base.getindex(b::LabelBasis, i) = tuple_at_ind(b, i)
    Base.getindex(b::LabelBasis, t::Tuple) = getpos(b, t)
    Base.getindex(b::LabelBasis, arr::AbstractArray) = [[b[i] for i in arr]...]

    ######################
    # Iterator Functions #
    ######################
    Base.start(::LabelBasis) = 1
    Base.done(b::LabelBasis, state) = length(b) == state-1
    Base.next(b::LabelBasis, state) = b[state], state+1
    Base.endof(b::LabelBasis) = length(b)
    Base.last(b::LabelBasis) = b[length(b)]
    Base.first(b::LabelBasis) = b[1]
    Base.collect(b::LabelBasis) = b[1:end]

    ##########################
    # Mathematical Functions #
    ##########################
    tensor{S}(a::LabelBasis{S}, b::LabelBasis{S}) = LabelBasis(tuple(a.ranges..., b.ranges...), S)

    ######################
    # Printing Functions #
    ######################
    Base.repr(b::LabelBasis) = "$(typeof(b))$(ranges(b))"
    Base.show(io::IO, b::LabelBasis) = print(io, repr(b))

export LabelBasis
