import Base: 
	getindex,
	length,
	convert,
	repr,
	show

##############
# StateLabel #
##############
	# A foundational label structure for Dirac states/operators 
	# Its flatten methods allow it to easily support tensor
	# product structures.
	# Note that objects of this type are not 
	# in and of themselves quantum objects.
	immutable StateLabel
		label::Tuple
		StateLabel(label::Tuple) = new(label)
		StateLabel(label...) = StateLabel(label)
	end

	StateLabel(s::StateLabel) = StateLabel(s.label)

	convert(::Type{StateLabel}, t::Tuple) = StateLabel(t)

#############
# Functions #
#############
	label(s::StateLabel) = s

	binarycombine(a::StateLabel, b::StateLabel) = StateLabel(tuple(a.label..., b.label...)) 
	combine(s::(StateLabel...)) = reduce(binarycombine, s)
	combine(s::StateLabel...) = combine(s)

	getindex(s::StateLabel, i) = getindex(s.label, i)
	length(s::StateLabel) = length(s.label)

	labelstr(label::String) = strip(match(r"\(.*?\)", label).match[2:end-1], ',')
	labelstr(s::StateLabel) = "$(labelstr(repr(s.label)))"
	repr(s::StateLabel) = "($(labelstr(s)))"
	show(io::IO, s::StateLabel) = print(io, repr(s))

export StateLabel,
	label,
	labelstr,
	combine