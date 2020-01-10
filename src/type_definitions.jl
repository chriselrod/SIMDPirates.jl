
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

@generated function vsetindex(v::AbstractSIMDVector{N,T}, x::Number, ::Val{I}) where {N,T,I}
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

@generated function vsetindex(v::AbstractSIMDVector{N,T}, x::Number, i::Int) where {N,T}
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

@inline Base.setindex(v::AbstractStructVec{N,T}, x::Number, i) where {N,T} = setindex(v, x, i)
@inline Base.setindex(v::AbstractStructVec{N,T}, x::Number, i::Integer) where {N,T} = setindex(v, x, Int(i))


@inline getvalindex(v::AbstractSIMDVector{N,T}, ::Val{I}) where {N,T,I} = extract_data(v)[I].value




@inline Vec{N,T}(v::Vararg{T,N}) where {T,N} = ntuple(n -> VE(v[n]), Val(N))

@inline Base.convert(::Type{SVec{W,T}}, v::Vec{W,T}) where {W,T} = SVec(v)
@generated function Base.convert(::Type{SVec{N,T}}, xs::NTuple{N,T}) where {N,T}
    svecq = Expr(:call, :SVec, Expr(:tuple, [:(VE(xs[$n])) for n ∈ 1:N]...))
    quote
        $(Expr(:meta,:inline))
        @inbounds $svecq
    end
end
@generated function Base.convert(::Type{SVec{N,T2}}, xs::NTuple{N,T1}) where {N,T1,T2}
    quote
        $(Expr(:meta,:inline))
        @inbounds SVec((Base.Cartesian.@ntuple $N n -> VE(convert(T2,xs[n]))))
    end
end
@generated function Base.convert(::Type{SVec{N,T2}}, xs::Vec{N,T1}) where {N,T1,T2}
    quote
        $(Expr(:meta,:inline))
        @inbounds SVec((Base.Cartesian.@ntuple $N n -> VE(convert(T2,xs[n].value))))
    end
end
@generated function Base.convert(::Type{SVec{N,T1}}, xs::SVec{N,T2}) where {N,T1,T2}
    quote
        $(Expr(:meta,:inline))
        @inbounds SVec((Base.Cartesian.@ntuple $N n -> VE{T1}(xs[n])))
    end
end
@inline Base.convert(::Type{SVec{N,T}}, x::SVec{N,T}) where {N,T} = x

@generated function vconvert(::Type{Vec{W,T1}}, v::Vec{W,T2}) where {W,T1 <: FloatingTypes, T2 <: Signed}
    typ1 = llvmtype(T1)
    typ2 = llvmtype(T2)
    vtyp1 = "<$W x $typ1>"
    vtyp2 = "<$W x $typ2>"
    instrs = """
%res = sitofp $vtyp2 %0 to $vtyp1
ret $vtyp1 %res
    """
    quote
        $(Expr(:meta, :inline))
        Base.llvmcall($instrs, Vec{$W,$T1}, Tuple{Vec{$W,$T2}}, v)
    end
end
@generated function vconvert(::Type{Vec{W,T1}}, v::Vec{W,T2}) where {W,T1 <: FloatingTypes, T2 <: Unsigned}
    typ1 = llvmtype(T1)
    typ2 = llvmtype(T2)
    vtyp1 = "<$W x $typ1>"
    vtyp2 = "<$W x $typ2>"
    instrs = """
%res = uitofp $vtyp2 %0 to $vtyp1
ret $vtyp1 %res
    """
    quote
        $(Expr(:meta, :inline))
        Base.llvmcall($instrs, Vec{$W,$T1}, Tuple{Vec{$W,$T2}}, v)
    end
end
@generated function vconvert(::Type{Vec{W,T1}}, v::Vec{W,T2}) where {W,T1 <: Signed, T2 <: FloatingTypes}
    typ1 = llvmtype(T1)
    typ2 = llvmtype(T2)
    vtyp1 = "<$W x $typ1>"
    vtyp2 = "<$W x $typ2>"
    instrs = """
%res = fptosi $vtyp2 %0 to $vtyp1
ret $vtyp1 %res
    """
    quote
        $(Expr(:meta, :inline))
        Base.llvmcall($instrs, Vec{$W,$T1}, Tuple{Vec{$W,$T2}}, v)
    end
end
@generated function vconvert(::Type{Vec{W,T1}}, v::Vec{W,T2}) where {W,T1 <: Unsigned, T2 <: FloatingTypes}
    typ1 = llvmtype(T1)
    typ2 = llvmtype(T2)
    vtyp1 = "<$W x $typ1>"
    vtyp2 = "<$W x $typ2>"
    instrs = """
%res = fptoui $vtyp2 %0 to $vtyp1
ret $vtyp1 %res
    """
    quote
        $(Expr(:meta, :inline))
        Base.llvmcall($instrs, Vec{$W,$T1}, Tuple{Vec{$W,$T2}}, v)
    end
end
@inline vconvert(::Type{Vec{W,T}}, v::SVec) where {W,T} = vconvert(Vec{W,T}, extract_data(v))
@inline vconvert(::Type{SVec{W,T}}, v) where {W,T} = SVec(vconvert(Vec{W,T}, extract_data(v)))
@inline vconvert(::Type{T}, v::T) where {T} = v
@inline vconvert(::Type{Vec{W,T}}, v::Vec{W,T}) where {W,T} = v
@inline vconvert(::Type{SVec{W,T}}, v::SVec{W,T}) where {W,T} = v
@inline Base.convert(::Type{SVec{W,T}}, v) where {W,T} = SVec(vconvert(Vec{W,T}, extract_data(v)))


@inline promote_vtype(::Type{T}, ::Type{T}) where {T} = T
@inline function promote_vtype(::Type{V1}, ::Type{V2}) where {W,T1,T2,V1<:AbstractSIMDVector{W,T1},V2<:AbstractSIMDVector{W,T2}}
    T = promote_type(T1, T2)
    Vec{W,T}
end
@inline function promote_vtype(::Type{SVec{W,T1}}, ::Type{V2}) where {W,T1,T2,V2<:AbstractSIMDVector{W,T2}}
    T = promote_type(T1, T2)
    SVec{W,T}
end
@inline function promote_vtype(::Type{V1}, ::Type{SVec{W,T2}}) where {W,T1,T2,V1<:AbstractSIMDVector{W,T1}}
    T = promote_type(T1, T2)
    SVec{W,T}
end
@inline function promote_vtype(::Type{SVec{W,T1}}, ::Type{SVec{W,T2}}) where {W,T1,T2}
    T = promote_type(T1, T2)
    SVec{W,T}
end

@generated function zeropad(v::Vec{W,T}) where {W,T}
    typ = llvmtype(T)
    W2 = W << 1
    sv = join((i for i ∈ 0:W2-1), ", i32 ")
    instr = """
    %res = shufflevector <$W x $typ> %0, <$W x $typ> zeroinitializer, <$W2 x i32> <i32 $(sv)>
    ret <$W2 x $typ> %res
"""
    :(Base.llvmcall($instr, Vec{$W2,$T}, Tuple{Vec{$W,$T}}, v))
end


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

@generated function vecbool_to_unsigned(vb::Vec{N,Bool}) where {N}
    Nout = VectorizationBase.nextpow2(max(N, 8))
    Utype = VectorizationBase.mask_type(Nout)
    btype = llvmtype(Bool)
    vbtyp = "<$N x $btype>"
    decls = String[]
    instrs = String["%maskb = trunc $vbtyp %0 to <$N x i1>"]
    u1name = Nout == N ? "%uout" : "%uabridged"
    push!(instrs, "$u1name = bitcast <$N x i1> %maskb to i$N")
    if Nout > N
        push!(instrs, "%uout = zext i$N %uabridged to i$Nout")
    end
    push!(instrs, "ret i$Nout %uout")
    quote
        $(Expr(:meta,:inline))
        Base.llvmcall($(join(instrs, "\n")), $Utype, Tuple{Vec{$N,Bool}}, vb)
    end
end

