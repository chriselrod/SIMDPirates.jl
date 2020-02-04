
const BoolTypes = Union{Bool}
const IntTypes = Union{Int8, Int16, Int32, Int64, Int128}
const UIntTypes = Union{UInt8, UInt16, UInt32, UInt64, UInt128}
const IntegerTypes = Union{IntTypes, UIntTypes, Ptr}
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

@inline Base.setindex(v::SVec{N,T}, x::Number, i) where {N,T} = setindex(v, x, i)
@inline Base.setindex(v::SVec{N,T}, x::Number, i::Integer) where {N,T} = setindex(v, x, Int(i))


@inline getvalindex(v::AbstractSIMDVector{N,T}, ::Val{I}) where {N,T,I} = extract_data(v)[I].value
@inline Base.convert(::Type{SVec{W,T}}, v::Vec{W,T}) where {W,T} = SVec(v)
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
@generated function vconvert(::Type{Vec{W,I1}}, v::Vec{W,I2}) where {W,I1<:Integer,I2<:Integer}
    typ1 = llvmtype(I1)
    typ2 = llvmtype(I2)
    b1 = 8sizeof(I1)
    b2 = 8sizeof(I2)
    op = if b1 < b2
        "trunc"
    elseif b1 > b2
        "zext"
    else
        return I1 === I2 ? :v : Expr(:block, Expr(:meta,:inline), Expr(:call, :vreinterpret, Vec{W,I1}, :v))
    end
    instrs = """
    %res = $op <$W x $typ2> %0 to <$W x $typ1>
    ret <$W x $typ1> %res
    """
    quote
        $(Expr(:meta, :inline))
        Base.llvmcall($instrs, Vec{$W,$I1}, Tuple{Vec{$W,$I2}}, v)
    end
end
@generated function vconvert(::Type{Vec{W,T1}}, v::Vec{W,T2}) where {W, T1 <: FloatingTypes, T2 <: FloatingTypes}
    typ1 = llvmtype(T1)
    typ2 = llvmtype(T2)
    b1 = 8sizeof(T1)
    b2 = 8sizeof(T2)
    op = if b1 < b2
        "fptrunc"
    elseif b1 > b2
        "fpext"
    end
    instrs = """
    %res = $op <$W x $typ2> %0 to <$W x $typ1>
    ret <$W x $typ1> %res
    """
    quote
        $(Expr(:meta, :inline))
        Base.llvmcall($instrs, Vec{$W,$T1}, Tuple{Vec{$W,$T2}}, v)
    end
end
@inline vconvert(::Type{Vec{W,T}}, v::Vec{W,T}) where {W,T<:Integer} = v # specific definition
@inline vconvert(::Type{Vec{W,T}}, v::SVec) where {W,T} = vconvert(Vec{W,T}, extract_data(v))
@inline vconvert(::Type{SVec{W,T}}, v) where {W,T} = SVec(vconvert(Vec{W,T}, extract_data(v)))
@inline vconvert(::Type{T}, v::T) where {T} = v
@inline vconvert(::Type{Vec{W,T}}, v::Vec{W,T}) where {W,T<:FloatingTypes} = v
@inline vconvert(::Type{Vec{W,T}}, v::Vec{W,T}) where {W,T} = v
@inline vconvert(::Type{SVec{W,T}}, v::SVec{W,T}) where {W,T} = v
@inline Base.convert(::Type{SVec{W,T}}, v) where {W,T} = SVec(vconvert(Vec{W,T}, extract_data(v)))
@inline VectorizationBase.SVec{W,T1}(v::SVec{W,T2}) where {W,T1<:FloatingTypes,T2<:IntegerTypes} = vconvert(SVec{W,T1}, v)
@inline VectorizationBase.SVec{W,T1}(v::SVec{W,T2}) where {W,T1<:IntegerTypes,T2<:FloatingTypes} = vconvert(SVec{W,T1}, v)
@inline vconvert(::Type{Vec{W,T1}}, s::T2) where {W, T1, T2 <: FloatingTypes} = vbroadcast(Vec{W,T1}, convert(T1, s))
@inline vconvert(::Type{Vec{W,T1}}, s::T2) where {W, T1 <: FloatingTypes, T2 <: Integer} = vbroadcast(Vec{W,T1}, convert(T1, s))
@inline vconvert(::Type{Vec{W,T1}}, s::T2) where {W, T1 <: Integer, T2 <: Integer} = vbroadcast(Vec{W,T1}, Base.unsafe_trunc(T1, s))
@inline vconvert(::Type{Vec{W,T1}}, s::T1) where {W, T1 <: Integer} = vbroadcast(Vec{W,T1}, s)

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
# @inline promote_vtype(::Type{V}, ::Type{T}) where {W,T,V <: AbstractSIMDVector{W,T}} = V
@inline function promote_vtype(::Type{Vec{W,T1}}, ::Type{T2}) where {W,T1,T2<:Number}
    T = promote_type(T1,T2)
    Vec{W,T}
end
@inline function promote_vtype(::Type{SVec{W,T1}}, ::Type{T2}) where {W,T1,T2<:Number}
    T = promote_type(T1,T2)
    SVec{W,T}
end
@inline function promote_vtype(::Type{T2}, ::Type{Vec{W,T1}}) where {W,T1,T2<:Number}
    T = promote_type(T1,T2)
    Vec{W,T}
end
@inline function promote_vtype(::Type{T2}, ::Type{SVec{W,T1}}) where {W,T1,T2<:Number}
    T = promote_type(T1,T2)
    SVec{W,T}
end
@inline promote_vtype(::Type{T1}, ::Type{T2}, ::Type{T3}) where {T1,T2,T3} = promote_vtype(promote_vtype(T1, T2), T3)
@inline promote_vtype(::Type{T}, ::Type{T}, ::Type{T}) where {T} = T

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

