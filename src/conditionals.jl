
# Conditionals

for op ∈ (:(==), :(!=), :(<), :(<=), :(>), :(>=) )
    rename = VECTOR_SYMBOLS[op]
    erename = Symbol(:e, rename)
    # rename_mask = Symbol(rename, :_mask)
    @eval begin
        # scalar versions handled in floating_point_arithmetic.jl
        # @inline $rename(s1::ScalarTypes, s2::ScalarTypes) = $op(s1,s2)
        @vectordef $rename function Base.$op(v1, v2) where {W,T}
            llvmwrap_bitmask(Val{$(QuoteNode(op))}(), extract_data(v1), extract_data(v2))
        end
        @vectordef $rename function Base.$op(v1, v2::T) where {W,T}
            $rename(extract_data(v1), v2)
        end
        @vectordef $rename function Base.$op(v1::T, v2) where {W,T}
            $rename(extract_data(v1), v2)
        end
        @evectordef $erename function Base.$op(v1, v2) where {W,T <: FloatingTypes}
            llvmwrap_bitmask_notfast(Val{$(QuoteNode(op))}(), extract_data(v1), extract_data(v2))
        end
    end
end
@inline visfinite(s::ScalarTypes) = isfinite(s)
@inline function visfinite(v1::Vec{W,T}) where {W,T<:FloatingTypes}
    U = uint_type(T)
    em = vbroadcast(Vec{W,U}, exponent_mask(T))
    iv = vreinterpret(Vec{W,U}, v1)
    evnot_equal(vand(iv, em), em)
end
@inline visfinite(v1::SVec{W}) where {W} = SVec{W}(visfinite(extract_data(v1)))
@inline Base.isfinite(v1::SVec{W}) where {W} = SVec{W}(visfinite(extract_data(v1)))

@inline visinf(s1::ScalarTypes) = isinf(s1)
@inline visinf(v1::Vec{W,T}) where {W,T<:FloatingTypes} = evisequal(vabs(v1), vbroadcast(Vec{W,T},Inf))
@inline visinf(v1::SVec{W}) where {W} = SVec{W}(visinf(extract_data(v1)))
@inline Base.isinf(v1::SVec{W}) where {W} = SVec{W}(visinf(extract_data(v1)))

@inline visnan(s1::ScalarTypes) = isnan(s1)
@inline visnan(v1::Vec{W,T}) where {W,T<:FloatingTypes} = evnot_equal(v1, v1)
@inline visnan(v1::SVec{W}) where {W} = SVec{W}(visnan(extract_data(v1)))
@inline Base.isnan(v1::SVec{W}) where {W} = SVec{W}(visnan(extract_data(v1)))

@inline vissubnormal(s1::ScalarTypes) = issubnormal(s1)
@inline function vissubnormal(v1::Vec{W,T}) where {W,T<:FloatingTypes}
    U = uint_type(T)
    em = vbroadcast(Vec{W,U}, exponent_mask(T))
    sm = vbroadcast(Vec{W,U}, significand_mask(T))
    iv = vreinterpret(Vec{W,U}, v1)
    vand(evisequal(vand(iv, em), vbroadcast(Vec{W,U}, 0)), evnot_equal(vand(iv, sm), vbroadcast(Vec{W,U}, 0)))
end
@inline vissubnormal(v1::SVec{W}) where {W} = SVec{W}(vissubnormal(extract_data(v1)))
@inline Base.issubnormal(v1::SVec{W}) where {W} = SVec{W}(vissubnormal(extract_data(v1)))


@inline signbit(s1::ScalarTypes) = signbit(s1)
@inline function vsignbit(v1::Vec{W,T}) where {W,T<:FloatingTypes}
    U = uint_type(T)
    sm = vbroadcast(Vec{W,U}, sign_mask(T))
    iv = vreinterpret(Vec{W,U}, v1)
    vnot_equal(vand(iv, sm), vbroadcast(Vec{W,U}, 0))
end
@inline function vsignbit(v1::SVec{W,T}) where {W,T<:FloatingTypes}
    SVec{W}(vsignbit(extract_data(v1)))
end
@inline function Base.signbit(v1::SVec{W,T}) where {W,T<:FloatingTypes}
    SVec{W}(vsignbit(extract_data(v1)))
end

@inline vifelse(c::Bool, x, y) = c ? x : y
# @inline vifelse(c::Bool, x, y) = ifelse(c, x, y)
@generated function vifelse(v1::Vec{W,Bool}, v2::Vec{W,T}, v3::Vec{W,T}) where {W,T}
    btyp = llvmtype(Bool)
    vbtyp = "<$W x $btyp>"
    abtyp = "[$W x $btyp]"
    typ = llvmtype(T)
    vtyp = "<$W x $typ>"
    atyp = "[$W x $typ]"
    decls = String[]
    instrs = String[]
    push!(instrs, "%cond = trunc $vbtyp %0 to <$W x i1>")
    if T <: FloatingTypes && Base.libllvm_version >= v"9"
        push!(instrs, "%res = select fast <$W x i1> %cond, $vtyp %1, $vtyp %2")
    else
        push!(instrs, "%res = select <$W x i1> %cond, $vtyp %1, $vtyp %2")
    end
    push!(instrs, "ret $vtyp %res")
    quote
        $(Expr(:meta, :inline))
        Base.llvmcall($((join(decls, "\n"), join(instrs, "\n"))),
            Vec{W,T},
            Tuple{Vec{W,Bool}, Vec{W,T}, Vec{W,T}},
            v1, v2, v3)
    end
end
@inline function vifelse(
    v1::AbstractSIMDVector{W,Bool}, v2::AbstractSIMDVector{W,T}, v3::AbstractSIMDVector{W,T}
) where {W,T}
    SVec(vifelse(extract_data(v1), extract_data(v2), extract_data(v3)))
end
# @inline function vifelse(v1::VecOrProd{W,Bool},
#         v2::VecOrProd{W,T}, v3::VecOrProd{W,T}) where {W,T}
#     vifelse(extract_data(v1), extract_data(v2), extract_data(v3))
# end
# @inline function Base.ifelse(v1::SVec{W,Bool},
#         v2::SVec{W,T}, v3::SVec{W,T}) where {W,T}
#     SVec(vifelse(extract_data(v1), extract_data(v2), extract_data(v3)))
# end

@generated function vifelse(mask::U, v2::Vec{W,T}, v3::Vec{W,T}) where {W,T,U<:Unsigned}
    @assert 8sizeof(U) >= W
    btyp = llvmtype(Bool)
    vbtyp = "<$W x $btyp>"
    abtyp = "[$W x $btyp]"
    typ = llvmtype(T)
    vtyp = "<$W x $typ>"
    atyp = "[$W x $typ]"
    decls = String[]
    instrs = String[]

    mtyp_input = llvmtype(U)
    mtyp_trunc = "i$W"

    if mtyp_input == mtyp_trunc
        push!(instrs, "%cond = bitcast $mtyp_input %0 to <$W x i1>")
    else
        push!(instrs, "%condtrunc = trunc $mtyp_input %0 to $mtyp_trunc")
        push!(instrs, "%cond = bitcast $mtyp_trunc %condtrunc to <$W x i1>")
    end
    # push!(instrs, "%cond = trunc $vbtyp %0 to <$W x i1>")
    if T <: FloatingTypes && Base.libllvm_version >= v"9"
        push!(instrs, "%res = select fast <$W x i1> %cond, $vtyp %1, $vtyp %2")
    else
        push!(instrs, "%res = select <$W x i1> %cond, $vtyp %1, $vtyp %2")
    end
    push!(instrs, "ret $vtyp %res")
    quote
        $(Expr(:meta, :inline))
        Base.llvmcall($((join(decls, "\n"), join(instrs, "\n"))),
            Vec{$W,$T},
            Tuple{$U, Vec{$W,$T}, Vec{$W,$T}},
            mask, v2, v3)
    end
end

@vpromote vifelse 3

@inline vifelse(U::Unsigned, v2::AbstractSIMDVector, v3::AbstractSIMDVector) = SVec(vifelse(U, extract_data(v2), extract_data(v3)))
@inline vifelse(U::Unsigned, v2::Vec{W,T}, s::Union{T,Int}) where {W,T} = vifelse(U, v2, vbroadcast(Vec{W,T}, s))
@inline vifelse(U::Unsigned, v2::AbstractSIMDVector{W,T}, s::Union{T,Int}) where {W,T} = SVec(vifelse(U, extract_data(v2), vbroadcast(Vec{W,T}, s)))
@inline vifelse(U::Unsigned, s::Union{T,Int}, v2::Vec{W,T}) where {W,T} = vifelse(U, vbroadcast(Vec{W,T}, s), v2)
@inline vifelse(U::Unsigned, s::Union{T,Int}, v2::AbstractSIMDVector{W,T}) where {W,T} = SVec(vifelse(U, vbroadcast(Vec{W,T}, s), extract_data(v2)))

@inline function vifelse(f::F, m::AbstractMask{W}, vargs::Vararg{<:Any,N}) where {F<:Function,W,N}
    vifelse(m, f(vargs...), @inbounds(vargs[N]))
end

@vectordef visodd function Base.isodd(v) where {W,T<:Integer}
    visequal(vand(v, one(T)), one(T))
end

# @inline vand(U::Unsigned, vb::Vec{W,Bool}) where {W} = U & vecbool_to_unsigned(vb)
# @inline vand(vb::Vec{W,Bool}, U::Unsigned) where {W} = U & vecbool_to_unsigned(vb)
