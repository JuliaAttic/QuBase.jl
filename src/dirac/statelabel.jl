import Base: getindex,
	length,
	convert,
	repr,
	show,
	start,
	done,
	next,
	endof,
	last,
	first,
	collect,
	map

##############
# StateLabel #
##############
	# The `StateLabel` type provides a foundational label structure 
	# for Dirac states/operators. The `combine` function allows it 
	# to easily support tensor product structures.
	# Note that objects of this type are not 
	# in and of themselves quantum objects.
	immutable StateLabel
		label::Tuple
	end
	StateLabel(label...) = StateLabel(label)
	StateLabel(s::StateLabel) = StateLabel(gettuple(s))
	StateLabel(s::AbstractState) = label(s)

	convert(::Type{StateLabel}, t::Tuple) = StateLabel(t)

	###############################
	# Accessor/Property Functions #
	###############################
	label(s::StateLabel) = s
	gettuple(s::StateLabel) = s.label
	getindex(s::StateLabel, i) = getindex(s.label, i)
	length(s::StateLabel) = length(s.label)

	#####################
	# Joining Functions #
	#####################
	binarycombine(a::StateLabel, b::StateLabel) = StateLabel(tuple(gettuple(a)..., gettuple(b)...)) 
	combine(s::(StateLabel...)) = reduce(binarycombine, s)
	combine(s::StateLabel...) = combine(s)

	separate(s::StateLabel) = map(StateLabel, gettuple(s))

	######################
	# Printing Functions #
	######################
	labelstr(label::Tuple) = strip(repr(label)[2:end-1], ',')
	labelstr(s::StateLabel) = "$(labelstr(s.label))"
	repr(s::StateLabel) = "StateLabel($(labelstr(s)))"
	show(io::IO, s::StateLabel) = print(io, repr(s))

	#################################
	# Iterator/Array-like Functions #
	#################################
	start(::StateLabel) = 1
	done(s::StateLabel, state) = length(s) == state-1
	next(s::StateLabel, state) = s[state], state+1
	endof(s::StateLabel) = length(s)
	last(s::StateLabel) = s[length(s)]
	first(s::StateLabel) = s[1]
	collect(s::StateLabel) = s[1:end]
	map(f::Union(Function,DataType), s::StateLabel) = StateLabel(map(f, gettuple(s)))
	map(f, s::StateLabel) = StateLabel(map(f, gettuple(s)))

	# using tuple() here allows arrays like [1, ('a', 'b'), :c]
	# to transform such that each element in the array acts as a 
	# single particle label (without it, the Tuple element would 
	# act like a multi-particle label)
	labelvec(A::AbstractArray) = [StateLabel(tuple(x)) for x in A]

export StateLabel,
	label,
	labelvec,
	gettuple,
	combine,
	separate