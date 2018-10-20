
const VE{T} = Core.VecElement{T}
const Vec{N,T} = NTuple{N,VE{T}}


const BoolTypes = Union{Bool}
const IntTypes = Union{Int8, Int16, Int32, Int64, Int128}
const UIntTypes = Union{UInt8, UInt16, UInt32, UInt64, UInt128}
const IntegerTypes = Union{BoolTypes, IntTypes, UIntTypes, Ptr}
const FloatingTypes = Union{Float16, Float32, Float64}
const ScalarTypes = Union{IntegerTypes, FloatingTypes}

const VE = Base.VecElement

type_length(::Type{Vec{N,T}}) where {N,T} = N
type_size(::Type{Vec{N,T}}) where {N,T} = (N,)
type_size(::Type{Vec{N,T}}, n::Integer) where {N,T} = (N,)[n]


# Element-wise access

export setindex
@generated function setindex(v::Vec{N,T}, x::Number, ::Type{Val{I}}) where {N,T,I}
    @assert isa(I, Integer)
    1 <= I <= N || throw(BoundsError())
    typ = llvmtype(T)
    ityp = llvmtype(Int)
    vtyp = "<$N x $typ>"
    decls = []
    instrs = []
    push!(instrs, "%res = insertelement $vtyp %0, $typ %1, $ityp $(I-1)")
    push!(instrs, "ret $vtyp %res")
    quote
        $(Expr(:meta, :inline))
        Vec{N,T}(Base.llvmcall($((join(decls, "\n"), join(instrs, "\n"))),
            NTuple{N,VE{T}}, Tuple{NTuple{N,VE{T}}, T}, v.elts, T(x)))
    end
end

@generated function setindex(v::Vec{N,T}, x::Number, i::Int) where {N,T}
    typ = llvmtype(T)
    ityp = llvmtype(Int)
    vtyp = "<$N x $typ>"
    decls = []
    instrs = []
    push!(instrs, "%res = insertelement $vtyp %0, $typ %2, $ityp %1")
    push!(instrs, "ret $vtyp %res")
    quote
        $(Expr(:meta, :inline))
        @boundscheck 1 <= i <= N || throw(BoundsError())
        Vec{N,T}(Base.llvmcall($((join(decls, "\n"), join(instrs, "\n"))),
            NTuple{N,VE{T}}, Tuple{NTuple{N,VE{T}}, Int, T},
            v.elts, i-1, T(x)))
    end
end

setindex(v::Vec{N,T}, x::Number, i) where {N,T} = setindex(v, Int(i), x)
getvalindex(v::Vec{N,T}, ::Type{Val{I}}) where {N,T,I} = v[I].value

@inline vbroadcast(::Type{Vec{N,T}}, s::S) where {N,T,S<:ScalarTypes} = ntuple(i -> VE{T}(s), Val(N))
@inline pirate_convert(::Type{Vec{N,T}}, xs::NTuple{N,T}) where {N,T<:ScalarTypes} = ntuple(i -> VE(xs[i]), Val(N))


# Floating point formats

int_type(::Type{Float16}) = Int16
int_type(::Type{Float32}) = Int32
int_type(::Type{Float64}) = Int64
# int_type(::Type{Float128}) = Int128
# int_type(::Type{Float256}) = Int256

uint_type(::Type{Float16}) = UInt16
uint_type(::Type{Float32}) = UInt32
uint_type(::Type{Float64}) = UInt64
# uint_type(::Type{Float128}) = UInt128
# uint_type(::Type{Float256}) = UInt256

significand_bits(::Type{Float16}) = 10
significand_bits(::Type{Float32}) = 23
significand_bits(::Type{Float64}) = 52
# significand_bits(::Type{Float128}) = 112
# significand_bits(::Type{Float256}) = 136

exponent_bits(::Type{T}) where {T<:FloatingTypes} =
    8*sizeof(T) - 1 - significand_bits(T)
sign_bits(::Type{T}) where {T<:FloatingTypes} = 1

significand_mask(::Type{T}) where {T<:FloatingTypes} =
    uint_type(T)(uint_type(T)(1) << significand_bits(T) - 1)
exponent_mask(::Type{T}) where {T<:FloatingTypes} =
    uint_type(T)(uint_type(T)(1) << exponent_bits(T) - 1) << significand_bits(T)
sign_mask(::Type{T}) where {T<:FloatingTypes} =
    uint_type(T)(1) << (significand_bits(T) + exponent_bits(T))

for T in (Float16, Float32, Float64)
    @assert sizeof(int_type(T)) == sizeof(T)
    @assert sizeof(uint_type(T)) == sizeof(T)
    @assert significand_bits(T) + exponent_bits(T) + sign_bits(T) == 8*sizeof(T)
    @assert significand_mask(T) | exponent_mask(T) | sign_mask(T) ==
        typemax(uint_type(T))
    @assert significand_mask(T) ⊻ exponent_mask(T) ⊻ sign_mask(T) ==
        typemax(uint_type(T))
end
