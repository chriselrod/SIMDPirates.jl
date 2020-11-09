
for T ∈ (:Float16,:Int16,:UInt16)
    @eval @inline sizeequivalentfloat(::Type{$T}) = Float16
    @eval @inline sizeequivalentint(::Type{$T}) = Int16
    @eval @inline sizeequivalentuint(::Type{$T}) = UInt16
    @eval @inline Base.$T(v::SVec{W}) where {W} = SVec(vconvert(Vec{W,$T}, extract_data(v)))
end
for T ∈ (:Float32,:Int32,:UInt32)
    @eval @inline sizeequivalentfloat(::Type{$T}) = Float32
    @eval @inline sizeequivalentint(::Type{$T}) = Int32
    @eval @inline sizeequivalentuint(::Type{$T}) = UInt32
    @eval @inline Base.$T(v::SVec{W}) where {W} = SVec(vconvert(Vec{W,$T}, extract_data(v)))
end
for T ∈ (:Float64,:Int64,:UInt64)
    @eval @inline sizeequivalentfloat(::Type{$T}) = Float64
    @eval @inline sizeequivalentint(::Type{$T}) = Int64
    @eval @inline sizeequivalentuint(::Type{$T}) = UInt64
    @eval @inline Base.$T(v::SVec{W}) where {W} = SVec(vconvert(Vec{W,$T}, extract_data(v)))
end
@inline sizeequivalentfloat(::Type{T}, x::T) where {T<:FloatingTypes} = x
@inline sizeequivalentint(::Type{T}, x::T) where {T<:Signed} = x
@inline sizeequivalentuint(::Type{T}, x::T) where {T<:Unsigned} = x
@inline function sizeequivalentfloat(::Type{T}, x) where {T}
    convert(sizeequivalentfloat(T), x)
end
@inline function sizeequivalentint(::Type{T}, x) where {T}
    convert(sizeequivalentint(T), x)
end
@inline function sizeequivalentint(::Type{T}, x::Integer) where {T}
    x % sizeequivalentint(T)
end
@inline function sizeequivalentuint(::Type{T}, x) where {T}
    convert(sizeequivalentint(T), x)
end
@inline function sizeequivalentuint(::Type{T}, x::Integer) where {T}
    x % sizeequivalentuint(T)
end

# Element-wise access

@generated function vsetindex(v::_Vec{_W,T}, x::Number, ::Val{I}) where {_W,T,I}
    W = _W + 1
    @assert isa(I, Integer)
    1 <= I <= W || throw(BoundsError())
    typ = llvmtype(T)
    ityp = llvmtype(Int)
    vtyp = "<$W x $typ>"
    # decls = String[]
    instrs = String[]
    push!(instrs, "%res = insertelement $vtyp %0, $typ %1, $ityp $(I-1)")
    push!(instrs, "ret $vtyp %res")
    quote
        $(Expr(:meta, :inline))
        llvmcall($(join(instrs, "\n")), Vec{$W,$T}, Tuple{Vec{$W,$T}, $T}, extract_data(v), $T(x))
    end
end

@generated function vsetindex(v::_Vec{_W,T}, x::Number, i::Int) where {_W,T}
    W = _W + 1
    typ = llvmtype(T)
    ityp = llvmtype(Int)
    vtyp = "<$W x $typ>"
    # decls = String[]
    instrs = String[]
    push!(instrs, "%res = insertelement $vtyp %0, $typ %2, $ityp %1")
    push!(instrs, "ret $vtyp %res")
    quote
        $(Expr(:meta, :inline))
        @boundscheck 1 <= i <= $W || throw(BoundsError())
        llvmcall($(join(instrs, "\n")), Vec{$W,$T}, Tuple{Vec{$W,$T}, Int, $T}, extract_data(v), vsub(i, 1), $T(x))
    end
end


Base.@propagate_inbounds Base.setindex(v::AbstractStructVec{W,T}, x::Number, i) where {W,T} = @inbounds SVec{W}(vsetindex(extract_data(v), x, i % Int))
@inline Base.setindex(v::AbstractStructVec{W,T}, x::Number, ::Val{I}) where {W,T,I} = SVec{W}(vsetindex(estract_data(v), x, I))



@inline getvalindex(v::AbstractSIMDVector, ::Val{I}) where {I} = extract_data(v)[I].value
@inline Base.convert(::Type{SVec{W,T}}, v::Vec{W,T}) where {W,T} = SVec(v)
@inline Base.convert(::Type{SVec{W,T}}, x::SVec{W,T}) where {W,T} = x

@inline VectorizationBase.vbroadcast(::Val{W}, v::AbstractSIMDVector{W}) where {W} = v
@inline vconvert(::Type{_Vec{_W,T}}, s::T) where {_W, T <: Unsigned} = vbroadcast(_Vec{_W,T}, s)
@generated function vconvert(::Type{_Vec{_W,T1}}, v::_Vec{_W,T2}) where {_W,T1 <: FloatingTypes, T2 <: Signed}
    W = _W + 1
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
        llvmcall($instrs, Vec{$W,$T1}, Tuple{Vec{$W,$T2}}, v)
    end
end
@generated function vconvert(::Type{_Vec{_W,T1}}, v::_Vec{_W,T2}) where {_W,T1 <: FloatingTypes, T2 <: Unsigned}
    W = _W + 1
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
        llvmcall($instrs, Vec{$W,$T1}, Tuple{Vec{$W,$T2}}, v)
    end
end
@generated function vconvert(::Type{_Vec{_W,T1}}, v::_Vec{_W,T2}) where {_W,T1 <: Signed, T2 <: FloatingTypes}
    W = _W + 1
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
        llvmcall($instrs, Vec{$W,$T1}, Tuple{Vec{$W,$T2}}, v)
    end
end
@generated function vconvert(::Type{_Vec{_W,T1}}, v::_Vec{_W,T2}) where {_W,T1 <: Unsigned, T2 <: FloatingTypes}
    W = _W + 1
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
        llvmcall($instrs, Vec{$W,$T1}, Tuple{Vec{$W,$T2}}, v)
    end
end
@generated function vconvert(::Type{_Vec{_W,I1}}, v::_Vec{_W,I2}) where {_W,I1<:Integer,I2<:Integer}
    W = _W + 1
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
        llvmcall($instrs, Vec{$W,$I1}, Tuple{Vec{$W,$I2}}, v)
    end
end
@generated function vconvert(::Type{_Vec{_W,T1}}, v::_Vec{_W,T2}) where {_W, T1 <: FloatingTypes, T2 <: FloatingTypes}
    W = _W + 1
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
        llvmcall($instrs, Vec{$W,$T1}, Tuple{Vec{$W,$T2}}, v)
    end
end
@inline vconvert(::Type{_Vec{_W,T}}, v::_Vec{_W,T}) where {_W,T<:Integer} = v # specific definition
@inline vconvert(::Type{_Vec{_W,T}}, v::SVec) where {_W,T} = vconvert(_Vec{_W,T}, extract_data(v))
@inline vconvert(::Type{SVec{W,T}}, v) where {W,T} = SVec(vconvert(Vec{W,T}, extract_data(v)))
@inline vconvert(::Type{T}, v::T) where {T} = v
@inline vconvert(::Type{_Vec{_W,T}}, v::_Vec{_W,T}) where {_W,T<:FloatingTypes} = v
@inline vconvert(::Type{_Vec{_W,T}}, v::_Vec{_W,T}) where {_W,T} = v
@inline vconvert(::Type{SVec{W,T}}, v::SVec{W,T}) where {W,T} = v
@inline Base.convert(::Type{SVec{W,T}}, v) where {W,T} = SVec(vconvert(Vec{W,T}, extract_data(v)))
@inline VectorizationBase.SVec{W,T1}(v::SVec{W,T2}) where {W,T1<:FloatingTypes,T2<:IntegerTypes} = vconvert(SVec{W,T1}, v)
@inline VectorizationBase.SVec{W,T1}(v::SVec{W,T2}) where {W,T1<:IntegerTypes,T2<:FloatingTypes} = vconvert(SVec{W,T1}, v)
@inline vconvert(::Type{_Vec{_W,T1}}, s::T2) where {_W, T1, T2 <: FloatingTypes} = vbroadcast(_Vec{_W,T1}, convert(T1, s))
@inline vconvert(::Type{_Vec{_W,T1}}, s::T2) where {_W, T1 <: FloatingTypes, T2 <: Integer} = vbroadcast(_Vec{_W,T1}, convert(T1, s))
@inline vconvert(::Type{_Vec{_W,T1}}, s::T2) where {_W, T1 <: Integer, T2 <: Integer} = vbroadcast(_Vec{_W,T1}, s % T1)
@inline vconvert(::Type{_Vec{_W,T1}}, s::T1) where {_W, T1 <: Integer} = vbroadcast(_Vec{_W,T1}, s)
@inline vconvert(::Type{T}, s::Number) where {T <: Number} = convert(T, s)
@inline vconvert(::Type{T}, s::Integer) where {T <: Integer} = s % T
@inline vconvert(::Type{<:Vec{W}}, m::Mask{W}) where {W} = m
@inline vconvert(::Type{<:Vec{W}}, u::Unsigned) where {W} = Mask{W}(u)
@inline vconvert(::Type{SVec{W,T}}, m::Mask{W}) where {W,T} = m
@inline vconvert(::Type{SVec{W,T1}}, s::T1) where {W, T1 <: Unsigned} = vbroadcast(SVec{W,T1}, s)
@inline vconvert(::Type{SVec{W,T}}, u::Unsigned) where {W,T} = Mask{W}(u)

@inline promote_vtype(::Type{T}, ::Type{T}) where {T} = T
@inline promote_vtype(::Type{Mask{W,U}}, ::Type{V}) where {W, U, V <: AbstractSIMDVector{W}} = V
@inline promote_vtype(::Type{V}, ::Type{Mask{W,U}}) where {W, U, V <: AbstractSIMDVector{W}} = V
@inline promote_vtype(::Type{Mask{W,U}}, ::Type{SVec{W,T}}) where {W, U, T} = SVec{W,T}
@inline promote_vtype(::Type{SVec{W,T}}, ::Type{Mask{W,U}}) where {W, U, T} = SVec{W,T}
@inline promote_vtype(::Type{Mask{W,U}}, ::Type{Mask{W,U}}) where {W, U} = Mask{W,U}
@inline promote_vtype(::Type{Mask{W,U}}, ::Type{T}) where {W, U, T <: Number} = SVec{W,T}
@inline promote_vtype(::Type{T}, ::Type{Mask{W,U}}) where {W, U, T <: Number} = SVec{W,T}
@inline promote_vtype(::Type{T1}, ::Type{T2}) where {T1 <: Number, T2 <: Number} = promote_type(T1, T2)
@inline promote_vtype(::Type{Vec{W,T}}, ::Type{SVec{W,T}}) where {W,T} = SVec{W,T}

@inline Base.unsafe_trunc(::Type{I}, v::SVec{W,T}) where {W, I <: Signed, T <: FloatingTypes} = SVec(vconvert(Vec{W,I}, extract_data(v)))

@inline function promote_vtype(::Type{V1}, ::Type{V2}) where {W,T1,T2,V1<:AbstractStructVec{W,T1},V2<:AbstractStructVec{W,T2}}
    T = promote_type(T1, T2)
    SVec{W,T}
end
@inline function promote_vtype(::Type{_Vec{W,T1}}, ::Type{_Vec{W,T2}}) where {W,T1,T2}
    T = promote_type(T1, T2)
    _Vec{W,T}
end
@inline promote_type_vs(::Type{V}, ::Type{S}) where {V,S} = V # fallback
@inline promote_type_vs(::Type{V}, ::Type{S}) where {V<:Signed,S<:Signed} = V
@inline promote_type_vs(::Type{V}, ::Type{S}) where {V<:Unsigned,S<:Signed} = V
@inline promote_type_vs(::Type{V}, ::Type{S}) where {V<:Signed,S<:Unsigned} = unsigned(V)
@inline promote_type_vs(::Type{V}, ::Type{S}) where {V<:Unsigned,S<:Unsigned} = V
@inline promote_type_vs(::Type{V}, ::Type{S}) where {V<:FloatingTypes,S<:FloatingTypes} = V
@inline promote_type_vs(::Type{V}, ::Type{S}) where {V<:Unsigned,S<:FloatingTypes} = Base.uinttype(S)
@inline promote_type_vs(::Type{V}, ::Type{S}) where {V<:Signed,S<:FloatingTypes} = signed(Base.uinttype(S))
@inline promote_type_vs(::Type{<:Integer}, ::Type{Ptr{T}}) where {T} = Ptr{T}
    
@inline function promote_vtype(::Type{_Vec{W,T1}}, ::Type{T2}) where {W,T1,T2<:Number}
    T = promote_type_vs(T1,T2)
    _Vec{W,T}
end
@inline function promote_vtype(::Type{SVec{W,T1}}, ::Type{T2}) where {W,T1,T2<:Number}
    T = promote_type_vs(T1,T2)
    SVec{W,T}
end
@inline function promote_vtype(::Type{T2}, ::Type{_Vec{W,T1}}) where {W,T1,T2<:Number}
    T = promote_type_vs(T1,T2)
    _Vec{W,T}
end
@inline function promote_vtype(::Type{T2}, ::Type{SVec{W,T1}}) where {W,T1,T2<:Number}
    T = promote_type_vs(T1,T2)
    SVec{W,T}
end
@inline promote_vtype(::Type{<:_MM{W}}, ::Type{T}) where {W,T<:Number} = SVec{W,T}
@inline promote_vtype(::Type{T}, ::Type{<:_MM{W}}) where {W,T<:Number} = SVec{W,T}
@inline promote_vtype(::Type{<:_MM{W}}, ::Type{V}) where {W,V<:Vec{W}} = V
@inline promote_vtype(::Type{V}, ::Type{<:_MM{W}}) where {W,V<:Vec{W}} = V
@inline promote_vtype(::Type{<:_MM{W}}, ::Type{SVec{W,T}}) where {W,T<:Number} = SVec{W,T}
@inline promote_vtype(::Type{SVec{W,T}}, ::Type{<:_MM{W}}) where {W,T<:Number} = SVec{W,T}
@inline promote_vtype(::Type{T1}, ::Type{T2}, ::Type{T3}) where {T1,T2,T3} = promote_vtype(promote_vtype(T1, T2), T3)
@inline promote_vtype(::Type{T}, ::Type{T}, ::Type{T}) where {T} = T
@inline promote_vtype(::Type{_Vec{W,T}}, ::Type{_Vec{W,T}}) where {W,T} = _Vec{W,T}
@inline promote_vtype(::Type{SVec{W,T}}, ::Type{SVec{W,T}}) where {W,T} = SVec{W,T}

@inline function vpromote(args::Vararg{<:Any,N}) where {N}
    V = promote_vtype(typeof.(args)...)
    vconvert.(V, args)
end

@inline Base.float(v::AbstractStructVec{W,<:FloatingTypes}) where {W} = v
@inline Base.float(v::AbstractStructVec{W,I}) where {W, I <: Union{Int64,UInt64}} = vconvert(SVec{W,Float64}, v)
@inline Base.float(v::AbstractStructVec{W,I}) where {W, I <: Union{Int32,UInt32}} = vconvert(SVec{W,Float32}, v)

@inline Base.unsigned(v::AbstractStructVec{W,Int8}) where {W} = vconvert(SVec{W,UInt8}, v)
@inline Base.unsigned(v::AbstractStructVec{W,Int16}) where {W} = vconvert(SVec{W,UInt16}, v)
@inline Base.unsigned(v::AbstractStructVec{W,Int32}) where {W} = vconvert(SVec{W,UInt32}, v)
@inline Base.unsigned(v::AbstractStructVec{W,Int64}) where {W} = vconvert(SVec{W,UInt64}, v)
@inline Base.unsigned(v::AbstractStructVec{W,Int128}) where {W} = vconvert(SVec{W,UInt128}, v)

@generated function zeropad(v::_Vec{_W,T}) where {_W,T}
    W = _W + 1
    typ = llvmtype(T)
    W2 = W << 1
    sv = join((i for i ∈ 0:W2-1), ", i32 ")
    instr = """
    %res = shufflevector <$W x $typ> %0, <$W x $typ> zeroinitializer, <$W2 x i32> <i32 $(sv)>
    ret <$W2 x $typ> %res
    """
    :(llvmcall($instr, Vec{$W2,$T}, Tuple{Vec{$W,$T}}, v))
end

@inline Base.trunc(::Type{I}, v::SVec{W,T}) where {W,I<:Integer,T<:FloatingTypes} = vconvert(SVec{W,I}, v)


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


# @generated function vecbool_to_unsigned(vb::Vec{N,Bool}) where {N}
#     Nout = VectorizationBase.nextpow2(max(N, 8))
#     Utype = VectorizationBase.mask_type(Nout)
#     btype = llvmtype(Bool)
#     vbtyp = "<$N x $btype>"
#     decls = String[]
#     instrs = String["%maskb = trunc $vbtyp %0 to <$N x i1>"]
#     u1name = Nout == N ? "%uout" : "%uabridged"
#     push!(instrs, "$u1name = bitcast <$N x i1> %maskb to i$N")
#     if Nout > N
#         push!(instrs, "%uout = zext i$N %uabridged to i$Nout")
#     end
#     push!(instrs, "ret i$Nout %uout")
#     quote
#         $(Expr(:meta,:inline))
#         llvmcall($(join(instrs, "\n")), $Utype, Tuple{Vec{$N,Bool}}, vb)
#     end
# end

