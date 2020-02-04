
# Integer arithmetic functions

for op ∈ (:(~), :(+), :(-))
    rename = VECTOR_SYMBOLS[op]
    @eval begin
        @inline $rename(s1::IntegerTypes) = $op($IntegerTypes)
        @vectordef $rename function Base.$op(v1) where {N,T<:IntegerTypes}
            llvmwrap(Val{$(QuoteNode(op))}(), extract_data(v1))
        end
    end
end
# @inline vnot(s1::Bool) = !s1
# @inline vnot(v1::Vec{N,Bool}) where {N} = vbitwise_not(v1)
# @inline vnot(v1::SVec{N,Bool}) where {N} = SVec(vbitwise_not(extract_data(v1)))
# @inline Base.:!(v1::SVec{N,Bool}) where {N} = SVec(vbitwise_not(extract_data(v1)))

@inline vabs(s1::IntTypes) = abs(s1)
@inline function vabs(v1::Vec{N,T}) where {N,T<:IntTypes}
    # s = -vbroadcast(Vec{N,T}, signbit(v1))
    s = vright_bitshift(v1, Val{8*sizeof(T)}())
    # Note: -v1 == ~v1 + 1
    vsub(vxor(s, v1), s)
end
@inline vabs(v1::SVec{N,T}) where {N,T<:IntTypes} = SVec(vabs(extract_data(v1)))
@inline Base.abs(v1::SVec{N,T}) where {N,T<:IntTypes} = SVec(vabs(extract_data(v1)))

@inline vabs(v1::AbstractSIMDVector{N,T}) where {N,T<:UIntTypes} = v1
@inline Base.abs(v1::SVec{N,T}) where {N,T<:UIntTypes} = v1

# TODO: Try T(v1>0) - T(v1<0)
#       use a shift for v1<0
#       evaluate v1>0 as -v1<0 ?

@inline vsign(s1::IntTypes) = sign(s1)
@inline vsign(v1::Vec{N,T}) where {N,T<:IntTypes} =
    vifelse(visequal(v1, vbroadcast(Vec{N,T},0)), vbroadcast(Vec{N,T},0),
        vifelse(vless(v1, vbroadcast(Vec{N,T},0)), vbroadcast(Vec{N,T},-1), vbroadcast(Vec{N,T},1)))
@inline vsign(v1::Vec{N,T}) where {N,T<:UIntTypes} =
    vifelse(visequal(v1, vbroadcast(Vec{N,T},0)), vbroadcast(Vec{N,T},0), vbroadcast(Vec{N,T},1))
    @inline vsign(v1::SVec) = SVec(vsign(extract_data(v1)))
    @inline Base.sign(v1::SVec) = SVec(vsign(extract_data(v1)))


@inline vsignbit(s1::IntegerTypes) = signbit(s1)
@inline vsignbit(v1::Vec{N,T}) where {N,T<:IntTypes} = vless(v1, vbroadcast(Vec{N,T}, 0))
@inline vsignbit(v1::SVec{N,T}) where {N,T<:IntTypes} = SVec(vsignbit(extract_data(v1)))
@inline Base.signbit(v1::SVec{N,T}) where {N,T<:IntTypes} = SVec(vsignbit(extract_data(v1)))


@inline vsignbit(v1::Vec{N,T}) where {N,T<:UIntTypes} = vbroadcast(Vec{N,Bool}, false)
@inline vsignbit(v1::SVec{N,T}) where {N,T<:UIntTypes} = vbroadcast(SVec{N,Bool}, false)
@inline Base.signbit(v1::SVec{N,T}) where {N,T<:UIntTypes} = vbroadcast(SVec{N,Bool}, false)

for op ∈ (:(&), :(|), :(⊻), :(+), :(-), :(*), :(÷), :(%) )
    rename = VECTOR_SYMBOLS[op]
    @eval begin
        @inline $rename(s1::IntegerTypes, s2::IntegerTypes) = $op(s1, s2)

        @vectordef $rename function Base.$op(v1, v2) where {N,T<:IntegerTypes}
            llvmwrap(Val{$(QuoteNode(op))}(), extract_data(v1), extract_data(v2))
        end

        # @inline $rename(v1::Vec{N,T}, v2::Vec{N,T}) where {N,T<:IntegerTypes} =
        #     llvmwrap(Val{$(QuoteNode(op))}(), v1, v2)
    end
end


@inline Base.:%(v::SVec{W,I}, ::Type{I}) where {W,I} = v
@inline Base.:%(v::SVec{W,I}, ::Type{Vec{W,I}}) where {W,I} = v
@inline Base.:%(v::SVec{W,I}, ::Type{SVec{W,I}}) where {W,I} = v
# @inline Base.:%(v::Vec{W,I}, ::Type{SVec{W,I}}) where {W,I} = SVec(v)
Base.:%(v::SVec{W,I1}, ::Type{I2}) where {W,I1,I2} =  SVec(vconvert(Vec{W,I2}, extract_data(v)))
@inline Base.:%(v::SVec{W,I1}, ::Type{SVec{W,I2}}) where {W,I1,I2} = v % I2
@inline Base.:%(v::SVec{W,I1}, ::Type{Vec{W,I2}}) where {W,I1,I2} = v % I2

@inline vrem(v1::V, ::Type{V}) where {V <: AbstractSIMDVector} = v1


@inline vcopysign(s1::IntegerTypes, s2::IntegerTypes) = copysign(s1, s2)
@inline vcopysign(v1::Vec{N,T}, v2::Vec{N,T}) where {N,T<:IntTypes} =
    vifelse(vsignbit(v2), -abs(v1), abs(v1))
# @inline function vcopysign(v1::SVec{N,T}, v2::SVec{N,T}) where {N,T}
    # SVec(vcopysign(extract_data(v1), extract_data(v2)))
# end
@inline function Base.copysign(v1::SVec{N,T}, v2::SVec{N,T}) where {N,T}
    SVec(vcopysign(extract_data(v1), extract_data(v2)))
end
@inline vcopysign(v1::AbstractSIMDVector{N,T}, v2::AbstractSIMDVector{N,T}) where {N,T<:UIntTypes} = v1
@inline Base.copysign(v1::SVec{N,T}, v2::SVec{N,T}) where {N,T<:UIntTypes} = v1

@inline vflipsign(s1::IntegerTypes, s2::IntegerTypes) = flipsign(s1, s2)
@inline vflipsign(v1::Vec{N,T}, v2::Vec{N,T}) where {N,T<:IntTypes} =
    vifelse(vsignbit(v2), -v1, v1)
@inline function vflipsign(v1::SVec{N,T}, v2::SVec{N,T}) where {N,T<:IntTypes}
    SVec(vflipsign(extract_data(v1), extract_data(V2)))
end
@inline function Base.flipsign(v1::SVec{N,T}, v2::SVec{N,T}) where {N,T<:IntTypes}
    SVec(vflipsign(extract_data(v1), extract_data(V2)))
end

@inline vflipsign(v1::AbstractSIMDVector{N,T}, v2::AbstractSIMDVector{N,T}) where {N,T<:UIntTypes} = v1
@inline Base.flipsign(v1::SVec{N,T}, v2::SVec{N,T}) where {N,T<:UIntTypes} = v1

@inline vmax(s1::IntegerTypes, s2::IntegerTypes) = vmax(s1, s2)
@inline vmax(v1::Vec{N,T}, v2::Vec{N,T}) where {N,T<:IntegerTypes} =
    vifelse(vless(v2, v1), v2, v1)
@inline vmin(v1::Vec{N,T}, v2::Vec{N,T}) where {N,T<:IntegerTypes} =
    vifelse(vless(v1, v2), v1, v2)

vmuladd(s1::ScalarTypes, s2::ScalarTypes, s3::ScalarTypes) = muladd(s1,s2,s3)
@inline function vmuladd(v1::Vec{N,T}, v2::Vec{N,T}, v3::Vec{N,T}) where {N,T<:IntegerTypes}
    vadd(vmul(v1,v2),v3)
end
@inline function vmuladd(v1::AbstractSIMDVector{N,T}, v2::AbstractSIMDVector{N,T},
        v3::AbstractSIMDVector{N,T}) where {N,T<:IntegerTypes}
    SVec(vadd(vmul(extract_data(v1),extract_data(v2)),extract_data(v3)))
end
@inline function Base.muladd(v1::SVec{N,T}, v2::SVec{N,T},
        v3::SVec{N,T}) where {N,T<:IntegerTypes}
    SVec(vadd(vmul(extract_data(v1),extract_data(v2)),extract_data(v3)))
end

vfma(s1::ScalarTypes, s2::ScalarTypes, s3::ScalarTypes) = fma(s1,s2,s3)
@inline function vfma(v1::Vec{N,T}, v2::Vec{N,T}, v3::Vec{N,T}) where {N,T<:IntegerTypes}
    vadd(vmul(v1,v2),v3)
end
@inline function vfma(v1::AbstractSIMDVector{N,T}, v2::AbstractSIMDVector{N,T},
        v3::AbstractSIMDVector{N,T}) where {N,T<:IntegerTypes}
    SVec(vadd(vmul(extract_data(v1),extract_data(v2)),extract_data(v3)))
end
@inline function Base.fma(v1::SVec{N,T}, v2::SVec{N,T},
        v3::SVec{N,T}) where {N,T<:IntegerTypes}
    SVec(vadd(vmul(extract_data(v1),extract_data(v2)),extract_data(v3)))
end

# TODO: Handle negative shift counts
#       use vifelse
#       ensure vifelse is efficient
for op ∈ (:(<<), :(>>), :(>>>))
    rename = VECTOR_SYMBOLS[op]
    @eval begin
        @inline $rename(s1::IntegerTypes, s2::IntegerTypes) = $op(s1, s2)

        @vectordef $rename function Base.$op(v1, ::Val{I}) where {N,T<:IntegerTypes,I}
            llvmwrapshift(Val{$(QuoteNode(op))}(), extract_data(v1), Val{I}())
        end
        @vectordef $rename function Base.$op(v1, x2::Unsigned) where {N,T<:IntegerTypes}
            llvmwrapshift(Val{$(QuoteNode(op))}(), extract_data(v1), x2)
        end
        @vectordef $rename function Base.$op(v1, x2::Int) where {N,T<:IntegerTypes}
            llvmwrapshift(Val{$(QuoteNode(op))}(), extract_data(v1), x2)
        end
        @vectordef $rename function Base.$op(v1, x2::Integer) where {N,T<:IntegerTypes}
            llvmwrapshift(Val{$(QuoteNode(op))}(), extract_data(v1), x2)
        end
        @vectordef $rename function Base.$op(v1, v2) where {N,T<:IntegerTypes}
            llvmwrapshift(Val{$(QuoteNode(op))}(), extract_data(v1), extract_data(v2))
        end
        @vectordef $rename function Base.$op(x1::T, v2) where {N,T<:IntegerTypes}
            $rename(vbroadcast(Vec{N,I}, x1), extract_data(v2))
        end
        @vectordef $rename function Base.$op(x1::I, v2) where {N,I<:Integer,T<:IntegerTypes}
            $rename(vbroadcast(Vec{N,I}, x1), vconvert(Vec{N,I},extract_data(v2)))
        end


        # @inline $rename(v1::Vec{N,T}, ::Val{I}) where {N,T<:IntegerTypes,I} =
            # llvmwrapshift(Val{$(QuoteNode(op))}(), v1, Val{I}())
        # @inline $rename(v1::Vec{N,T}, x2::Unsigned) where {N,T<:IntegerTypes} =
        #     llvmwrapshift(Val{$(QuoteNode(op))}(), v1, x2)
        # @inline $rename(v1::Vec{N,T}, x2::Int) where {N,T<:IntegerTypes} =
        #     llvmwrapshift(Val{$(QuoteNode(op))}(), v1, x2)
        # @inline $rename(v1::Vec{N,T}, x2::Integer) where {N,T<:IntegerTypes} =
        #     llvmwrapshift(Val{$(QuoteNode(op))}(), v1, x2)
        # @inline $rename(v1::Vec{N,T},
        #                  v2::Vec{N,U}) where {N,T<:IntegerTypes,U<:UIntTypes} =
        #     llvmwrapshift(Val{$(QuoteNode(op))}(), v1, v2)
        # @inline $rename(v1::Vec{N,T},
        #                  v2::Vec{N,U}) where {N,T<:IntegerTypes,U<:IntegerTypes} =
        #     llvmwrapshift(Val{$(QuoteNode(op))}(), v1, v2)
        # @inline $rename(x1::T, v2::Vec{N,T}) where {N,T<:IntegerTypes} =
        #     $rename(vbroadcast(Vec{N,T}, x1), v2)
    end
end


# Promote scalars of all IntegerTypes to vectors of IntegerTypes, leaving the
# vector type unchanged

for op ∈ (
        :(==), :(!=), :(<), :(<=), :(>), :(>=),
        :(&), :(|), :(⊻), :(+), :(-), :(*),
        :copysign, :div, :flipsign, :max, :min, :rem
    )
    rename = VECTOR_SYMBOLS[op]
    @eval begin
        @inline $rename(s1::ScalarTypes, s2::ScalarTypes) = $op(s1, s2)

        @vectordef $rename function Base.$op(s1::IntegerTypes, v2) where {N, T <: IntegerTypes}
            $rename(vbroadcast(Vec{N,T}, s1), extract_data(v2))
        end
        @vectordef $rename function Base.$op(v1, s2::IntegerTypes) where {N, T <: IntegerTypes}
            $rename(extract_data(v1), vbroadcast(Vec{N,T}, s2))
        end


        # @vectordef $rename function Base.$op(s1::Bool, v2) where {N, T <: Bool}
        #     $rename(vbroadcast(Vec{N,Bool},s1), extract_data(v2))
        # end
        # @vectordef $rename function Base.$op(v2, s1::Bool) where {N, T <: Bool}
        #     $rename(vbroadcast(Vec{N,Bool},s1), extract_data(v2))
        # end
        # @inline $rename(s1::Bool, v2::Vec{N,Bool}) where {N} =
        #     $op(vbroadcast(Vec{N,Bool},s1), v2)
        # @inline $rename(s1::IntegerTypes, v2::Vec{N,T}) where {N,T<:IntegerTypes} =
        #     $op(vbroadcast(Vec{N,T}, s1), v2)
        # @inline $rename(v1::Vec{N,T}, s2::IntegerTypes) where {N,T<:IntegerTypes} =
        #     $op(v1, vbroadcast(Vec{N,T}, s2))
        #
        # @inline $rename(s1::Bool, v2::SVec{N,Bool}) where {N} =
        #     SVec($op(Vec{N,Bool}(s1), extract_data(v2)))
        # @inline $rename(s1::IntegerTypes, v2::SVec{N,T}) where {N,T<:IntegerTypes} =
        #     SVec($op(vbroadcast(Vec{N,T}, s1), extract_data(v2)))
        # @inline $rename(v1::SVec{N,T}, s2::IntegerTypes) where {N,T<:IntegerTypes} =
        #     SVec($op(extract_data(v1), vbroadcast(Vec{N,T}, s2)))
    end
end
# @vectordef vifelse function Base.ifelse(c::Vec{N,Bool}, s1::T, v2) where {N,T<:IntegerTypes}
#     vifelse(c, vbroadcast(Vec{N,T}, s1), extract_data(v2))
# end
# @vectordef vifelse function Base.ifelse(c::Vec{N,Bool}, v1, s2::T) where {N,T<:IntegerTypes}
#     vifelse(c, extract_data(v1), vbroadcast(Vec{N,T}, s2))
# end
@inline function vifelse(c::Vec{N,Bool}, s1::T, v2::Vec{N,T}) where {N,T<:IntegerTypes}
    vifelse(c, vbroadcast(Vec{N,T}, s1), v2)
end
@inline function vifelse(c::Vec{N,Bool}, v1::Vec{N,T}, s2::T) where {N,T<:IntegerTypes}
    vifelse(c, v1, vbroadcast(Vec{N,T}, s2))
end
@inline function vifelse(c::AbstractSIMDVector{N,Bool}, s1::T, v2::AbstractSIMDVector{N,T}) where {N,T<:IntegerTypes}
    SVec(vifelse(extract_data(c), vbroadcast(Vec{N,T}, s1), extract_data(v2)))
end
@inline function vifelse(c::AbstractSIMDVector{N,Bool}, v1::AbstractSIMDVector{N,T}, s2::T) where {N,T<:IntegerTypes}
    SVec(vifelse(extract_data(c), extract_data(v1), vbroadcast(Vec{N,T}, s2)))
end

@inline evadd(v1::AbstractSIMDVector{W,T}, v2::AbstractSIMDVector{W,T}) where {W, T <: Integer} = vadd(v1, v2)
@inline evadd(v1::AbstractSIMDVector{W,T}, v2::T) where {W, T <: Integer} = vadd(v1, vbroadcast(Vec{W,T}(), v2))
@inline evadd(v1::T, v2::AbstractSIMDVector{W,T}) where {W, T <: Integer} = vadd(vbroadcast(Vec{W,T}(), v1), v2)
@inline evsub(v1::AbstractSIMDVector{W,T}, v2::AbstractSIMDVector{W,T}) where {W, T <: Integer} = vsub(v1, v2)
@inline evsub(v1::AbstractSIMDVector{W,T}, v2::T) where {W, T <: Integer} = vsub(v1, vbroadcast(Vec{W,T}(), v2))
@inline evsub(v1::T, v2::AbstractSIMDVector{W,T}) where {W, T <: Integer} = vsub(vbroadcast(Vec{W,T}(), v1), v2)
@inline evmul(v1::AbstractSIMDVector{W,T}, v2::AbstractSIMDVector{W,T}) where {W, T <: Integer} = vmul(v1, v2)
@inline evmul(v1::AbstractSIMDVector{W,T}, v2::T) where {W, T <: Integer} = vmul(v1, vbroadcast(Vec{W,T}(), v2))
@inline evmul(v1::T, v2::AbstractSIMDVector{W,T}) where {W, T <: Integer} = vmul(vbroadcast(Vec{W,T}(), v1), v2)


# @inline function vifelse(c::Vec{N,Bool}, s1::IntegerTypes, v2::Vec{N,T}) where {N,T<:IntegerTypes}
#     vifelse(c, vbroadcast(Vec{N,T}, s1), v2)
# end
# @inline function vifelse(c::Vec{N,Bool}, v1::Vec{N,T}, s2::IntegerTypes) where {N,T<:IntegerTypes}
#     vifelse(c, v1, vbroadcast(Vec{N,T}, s2))
# end

# for op ∈ (:muladd,)
#     rename = VECTOR_SYMBOLS[op]
#     @eval begin
#         @inline $rename(s1::IntegerTypes, v2::Vec{N,T},
#                 v3::Vec{N,T}) where {N,T<:IntegerTypes} =
#             $rename(vbroadcast(Vec{N,T}, s1), v2, v3)
#         @inline $rename(v1::Vec{N,T}, s2::IntegerTypes,
#                 v3::Vec{N,T}) where {N,T<:IntegerTypes} =
#             $rename(v1, vbroadcast(Vec{N,T}, s2), v3)
#         @inline $rename(s1::IntegerTypes, s2::IntegerTypes,
#                 v3::Vec{N,T}) where {N,T<:IntegerTypes} =
#             $rename(vbroadcast(Vec{N,T}, s1), vbroadcast(Vec{N,T}, s2), v3)
#         @inline $rename(v1::Vec{N,T}, v2::Vec{N,T},
#                 s3::IntegerTypes) where {N,T<:IntegerTypes} =
#             $rename(v1, v2, vbroadcast(Vec{N,T}, s3))
#         @inline $rename(s1::IntegerTypes, v2::Vec{N,T},
#                 s3::IntegerTypes) where {N,T<:IntegerTypes} =
#             $rename(vbroadcast(Vec{N,T}, s1), v2, vbroadcast(Vec{N,T}, s3))
#         @inline $rename(v1::Vec{N,T}, s2::IntegerTypes,
#                 s3::IntegerTypes) where {N,T<:IntegerTypes} =
#             $rename(v1, vbroadcast(Vec{N,T}, s2), vbroadcast(Vec{N,T}, s3))
#     end
# end

@generated function vpmaddwd(a::NTuple{W,Core.VecElement{Int16}}, b::NTuple{W,Core.VecElement{Int16}}) where {W}
    Wh = W >>> 1
    @assert 2Wh == W
    @assert (REGISTER_SIZE >> 1) ≥ W
    S = W * 16
    # decl = "@llvm.x86.avx512.pmaddw.d.512"
    instr = "@llvm.x86.avx512.pmaddw.d.$S"
    decl = "declare <$Wh x i32> $instr(<32 x i16>, <32 x i16>)"
    instrs = String[
        "%res = call <$Wh x i32> $instr(<$W x i16> %0, <$W x i16> %1)",
        "ret <$Wh x i32> %res"
    ]
    quote
        $(Expr(:meta,:inline))
        Base.llvmcall(
            $((decl,join(instrs,"\n"))),
            NTuple{$Wh,Core.VecElement{Int32}},
            Tuple{NTuple{$W,Core.VecElement{Int16}},NTuple{$W,Core.VecElement{Int16}}},
            a, b
        )
    end
end

@generated function vleading_zeros(v::NTuple{W,Core.VecElement{I}}) where {W, I <: Integer}
    typ = "i$(8sizeof(I))"
    vtyp = "<$W x $typ>"
    instr = "@llvm.ctlz.v$(W)$(typ)"
    decl = "declare $vtyp $instr($vtyp, i1)"
    instrs = String[
        "%res = call $vtyp $instr($vtyp %0, i1 1)"
        "ret $vtyp %res"
    ]
    quote
        $(Expr(:meta, :inline))
        Base.llvmcall(
            $((decl, join(instrs, "\n"))),
            Vec{$W,$I}, Tuple{Vec{$W,$I}}, v
        )
    end
end
@inline vleading_zeros(v::SVec{<:Any,<:Integer}) = SVec(vleading_zeros(extract_data(v)))
@inline Base.leading_zeros(v::SVec{<:Any,<:Integer}) = SVec(vleading_zeros(extract_data(v)))

