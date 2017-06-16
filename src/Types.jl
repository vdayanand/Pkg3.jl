module Types

using Base.Random: UUID
using Base.Intrinsics: ult_int, bswap_int

export SHA1, VersionSpec, @vs_str

## ordering of UUIDs ##

Base.isless(a::UUID, b::UUID) = a.value < b.value

## SHA1 ##

primitive type SHA1 160 end

function Base.convert(::Type{SHA1}, bytes::Vector{UInt8})
	length(bytes) == sizeof(SHA1) ||
		throw(ArgumentError("wrong number of bytes for SHA1 hash: $(length(bytes))"))
	return reinterpret(SHA1, bytes)[1]
end
Base.convert(::Type{Vector{UInt8}}, hash::SHA1) = reinterpret(UInt8, [hash])

Base.isless(a::SHA1, b::SHA1) = ult_int(bswap_int(a), bswap_int(b))

function Base.show(io::IO, hash::SHA1)
	print(io, "SHA1(")
	for octet in Vector{UInt8}(hash)
		print(io, hex(octet))
	end
	print(io, ")")
end

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
