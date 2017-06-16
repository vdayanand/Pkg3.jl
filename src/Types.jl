module Types

export VersionSpec, @vs_str

## VersionSpec ##

struct VersionSpec{n}
	spec::NTuple{n,Int}
	function VersionSpec{n}(spec::NTuple{n,Integer}) where n
		n <= 3 || throw(ArgumentError("VersionSpec: you can only specify major, minor and patch versions"))
		return new(spec)
	end
end
VersionSpec(spec::Integer...) = VersionSpec{length(spec)}(spec)

macro vs_str(s::String)
	isempty(s) || s == "*" ? VersionSpec() :
	VersionSpec(map(x->parse(Int, x), split(s, '.'))...)
end
function Base.show(io::IO, s::VersionSpec)
	print(io, "vs\"")
	join(io, s.spec, '.')
	print(io, '"')
end
Base.show(io::IO, ::VersionSpec{0}) = print(io, "vs\"*\"")
Base.getindex(s::VersionSpec, i::Int) = s.spec[i]

Base.in(v::VersionNumber, s::VersionSpec{0}) = true
Base.in(v::VersionNumber, s::VersionSpec{1}) = v.major == s.spec[1]
Base.in(v::VersionNumber, s::VersionSpec{2}) = v.major == s.spec[1] && v.minor == s.spec[2]
Base.in(v::VersionNumber, s::VersionSpec{3}) = v.major == s.spec[1] && v.minor == s.spec[2] && v.patch == s.spec[3]

end # module
