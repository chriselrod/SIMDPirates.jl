
const BoolTypes = Union{Bool}
const IntTypes = Union{Int8, Int16, Int32, Int64, Int128}
const UIntTypes = Union{UInt8, UInt16, UInt32, UInt64, UInt128}
const IntegerTypes = Union{BoolTypes, IntTypes, UIntTypes, Ptr}
const FloatingTypes = Union{Float16, Float32, Float64}
const ScalarTypes = Union{IntegerTypes, FloatingTypes}


type_length(::Type{Vec{N,T}}) where {N,T} = N
type_size(::Type{Vec{N,T}}) where {N,T} = (N,)
type_size(::Type{Vec{N,T}}, n::Integer) where {N,T} = (N,)[n]


# Element-wise access

export setindex
@generated function setindex(v::AbstractSIMDVector{N,T}, x::Number, ::Val{I}) where {N,T,I}
    @assert isa(I, Integer)
    1 <= I <= N || throw(BoundsError())
    typ = llvmtype(T)
    ityp = llvmtype(Int)
    vtyp = "<$N x $typ>"
    decls = String[]
    instrs = String[]
    push!(instrs, "%res = insertelement $vtyp %0, $typ %1, $ityp $(I-1)")
    push!(instrs, "ret $vtyp %res")
    quote
        $(Expr(:meta, :inline))
        Vec{N,T}(Base.llvmcall($((join(decls, "\n"), join(instrs, "\n"))),
            NTuple{N,VE{T}}, Tuple{NTuple{N,VE{T}}, T}, extract_data(v), T(x)))
    end
end

@generated function setindex(v::AbstractSIMDVector{N,T}, x::Number, i::Int) where {N,T}
    typ = llvmtype(T)
    ityp = llvmtype(Int)
    vtyp = "<$N x $typ>"
    decls = String[]
    instrs = String[]
    push!(instrs, "%res = insertelement $vtyp %0, $typ %2, $ityp %1")
    push!(instrs, "ret $vtyp %res")
    quote
        $(Expr(:meta, :inline))
        @boundscheck 1 <= i <= N || throw(BoundsError())
        Vec{N,T}(Base.llvmcall($((join(decls, "\n"), join(instrs, "\n"))),
            NTuple{N,VE{T}}, Tuple{NTuple{N,VE{T}}, Int, T},
            extract_data(v), i-1, T(x)))
    end
end

@inline setindex(v::AbstractSIMDVector{N,T}, x::Number, i) where {N,T} = setindex(v, Int(i), x)
@inline getvalindex(v::AbstractSIMDVector{N,T}, ::Val{I}) where {N,T,I} = extract_data(v)[I].value

@inline function vbroadcast(::Type{Vec{N,T}}, s::S) where {N,T,S<:ScalarTypes}
    @inbounds ntuple(i -> VE{T}(s), Val(N))
end
@inline function vbroadcast(::Type{Vec{N,T}}, s::Core.VecElement{T}) where {N,T}
    @inbounds ntuple(i -> s, Val(N))
end
@inline vbroadcast(::Type{Vec{N,T}}, s::Vec{N,T}) where {N,T} = s
@inline function vbroadcast(::Type{Vec{N,T}}, ptr::Ptr{T}) where {N,T}
    s = Core.VecElement{T}(VectorizationBase.load(ptr))
    @inbounds ntuple(_ -> s, Val(N))
end
@inline function vbroadcast(::Type{Vec{N,T}}, ptr::Ptr) where {N,T}
    s = Core.VecElement{T}(VectorizationBase.load(Base.unsafe_convert(Ptr{T},ptr)))
    @inbounds ntuple(_ -> s, Val(N))
end
@inline function svbroadcast(::Type{SVec{N,T}}, s::S) where {N,T,S<:ScalarTypes}
    @inbounds SVec(ntuple(i -> VE{T}(s), Val(N)))
end
@inline svbroadcast(::Type{SVec{N,T}}, s::SVec{N,T}) where {N,T} = s
@inline function vbroadcast(::Type{SVec{N,T}}, s::S) where {N,T,S<:ScalarTypes}
    @inbounds SVec(ntuple(i -> VE{T}(s), Val(N)))
end
@inline function vbroadcast(::Type{SVec{N,T}}, s::Core.VecElement{T}) where {N,T}
    @inbounds SVec(ntuple(i -> s, Val(N)))
end
@inline vbroadcast(::Type{SVec{N,T}}, s::Vec{N,T}) where {N,T} = SVec(s)
@inline vbroadcast(::Type{SVec{N,T}}, s::SVec{N,T}) where {N,T} = s
    
@inline vone(::Type{Vec{N,T}}) where {N,T} = ntuple(i -> Core.VecElement(one(T)), Val(N))
@inline vzero(::Type{Vec{N,T}}) where {N,T} = ntuple(i -> Core.VecElement(zero(T)), Val(N))
@inline vone(::Type{SVec{N,T}}) where {N,T} = SVec(ntuple(i -> Core.VecElement(one(T)), Val(N)))
@inline vzero(::Type{SVec{N,T}}) where {N,T} = SVec(ntuple(i -> Core.VecElement(zero(T)), Val(N)))
@inline vone(::Type{T}) where {T} = one(T)
@inline vzero(::Type{T}) where {T} = zero(T)

@inline Vec{N,T}(v::Vararg{T,N}) where {T,N} = ntuple(n -> VE(v[n]), Val(N))

@inline function pirate_convert(::Type{Vec{N,T}}, xs::NTuple{N,T}) where {N,T<:ScalarTypes}
    @inbounds ntuple(i -> VE(xs[i]), Val(N))
end
@inline function pirate_convert(::Type{Vec{N,T1}}, xs::Vec{N,T2}) where {N,T1<:ScalarTypes,T2<:ScalarTypes}
    @inbounds ntuple(i -> VE(T1(xs[i].value)), Val(N))
end
@generated function Base.convert(::Type{SVec{N,T}}, xs::NTuple{N,T}) where {N,T}
    quote
        $(Expr(:meta,:inline))
        @inbounds SVec((Base.Cartesian.@ntuple $N n -> VE(xs[n])))
    end
end
@generated function Base.convert(::Type{SVec{N,T1}}, xs::SVec{N,T2}) where {N,T1,T2}
    quote
        $(Expr(:meta,:inline))
        @inbounds SVec((Base.Cartesian.@ntuple $N n -> VE{T1}(xs[n])))
    end
end
@inline Base.convert(::Type{SVec{N,T}}, x::SVec{N,T}) where {N,T} = x

# @inline similar(s::T, ::Vec{N,T}, ::Vec{N,T}) = ntuple(i -> VE{T}(s), Val(N))
# @inline similar(s::T, ::Vec{N,T}, ::Vec{N,T}) = ntuple(i -> VE{T}(s), Val(N))

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
