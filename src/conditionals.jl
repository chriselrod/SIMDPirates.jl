
# Conditionals

for op ∈ (:(==), :(!=), :(<), :(<=), :(>), :(>=) )
    rename = VECTOR_SYMBOLS[op]
    # rename_mask = Symbol(rename, :_mask)
    rename_bool = Symbol(rename, :_boolvec)
    @eval begin
        # scalar versions handled in floating_point_arithmetic.jl
        # @inline $rename(s1::ScalarTypes, s2::ScalarTypes) = $op(s1,s2)
        @vectordef $rename function Base.$op(v1, v2) where {N,T}
            llvmwrap_bitmask(Val{$(QuoteNode(op))}(), extract_data(v1), extract_data(v2))
        end
        @vectordef $rename function Base.$op(v1, v2::T) where {N,T}
            $rename(extract_data(v1), v2)
        end
        @vectordef $rename function Base.$op(v1::T, v2) where {N,T}
            $rename(extract_data(v1), v2)
        end
        @evectordef $rename_bool function Base.$op(v1, v2) where {N,T}
            llvmwrap(Val{$(QuoteNode(op))}(), extract_data(v1), extract_data(v2), Bool)
        end
        @evectordef $rename_bool function Base.$op(v1, v2::T) where {N,T}
            $rename(extract_data(v1), v2)
        end
        @evectordef $rename_bool function Base.$op(v1::T, v2) where {N,T}
            $rename(extract_data(v1), v2)
        end
        # @evectordef $rename_mask function Base.$op(v1, v2) where {N,T}
        #     llvmwrap_bitmask(Val{$(QuoteNode(op))}(), extract_data(v1), extract_data(v2))
        # end
        # @evectordef $rename_mask function Base.$op(v1, v2::T) where {N,T}
        #     $rename(extract_data(v1), v2)
        # end
        # @evectordef $rename_mask function Base.$op(v1::T, v2) where {N,T}
        #     $rename(extract_data(v1), v2)
        # end
    end
end
@inline visfinite(s::ScalarTypes) = isfinite(s)
@inline function visfinite(v1::Vec{N,T}) where {N,T<:FloatingTypes}
    U = uint_type(T)
    em = vbroadcast(Vec{N,U}, exponent_mask(T))
    iv = vreinterpret(Vec{N,U}, v1)
    vnot_equal(vand(iv, em), em)
end
@inline visfinite(v1::SVec) = SVec(visfinite(extract_data(v1)))
@inline Base.isfinite(v1::SVec) = SVec(visfinite(extract_data(v1)))

@inline visfinite_boolvec(s::ScalarTypes) = isfinite(s)
@inline function visfinite_boolvec(v1::Vec{N,T}) where {N,T<:FloatingTypes}
    U = uint_type(T)
    em = vbroadcast(Vec{N,U}, exponent_mask(T))
    iv = vreinterpret(Vec{N,U}, v1)
    vnot_equal_boolvec(vand(iv, em), em)
end
@inline visfinite_boolvec(v1::SVec) = SVec(visfinite_boolvec(extract_data(v1)))


@inline visinf(s1::ScalarTypes) = isinf(s1)
@inline visinf(v1::Vec{N,T}) where {N,T<:FloatingTypes} = visequal(vabs(v1), vbroadcast(Vec{N,T},Inf))
@inline visinf(v1::SVec) = SVec(visinf(extract_data(v1)))
@inline Base.isinf(v1::SVec) = SVec(visinf(extract_data(v1)))

@inline visnan(s1::ScalarTypes) = isnan(s1)
@inline visnan(v1::Vec{N,T}) where {N,T<:FloatingTypes} = vnot_equal(v1, v1)
@inline visnan(v1::SVec) = SVec(visnan(extract_data(v1)))
@inline Base.isnan(v1::SVec) = SVec(visnan(extract_data(v1)))

@inline vissubnormal(s1::ScalarTypes) = issubnormal(s1)
@inline function vissubnormal(v1::Vec{N,T}) where {N,T<:FloatingTypes}
    U = uint_type(T)
    em = vbroadcast(Vec{N,U}, exponent_mask(T))
    sm = vbroadcast(Vec{N,U}, significand_mask(T))
    iv = vreinterpret(Vec{N,U}, v1)
    vand(visequal(vand(iv, em), vbroadcast(Vec{N,U}, 0)), vnot_equal(vand(iv, sm), vbroadcast(Vec{N,U}, 0)))
end
@inline vissubnormal(v1::SVec) = SVec(vissubnormal(extract_data(v1)))
@inline Base.issubnormal(v1::SVec) = SVec(vissubnormal(extract_data(v1)))


@inline signbit(s1::ScalarTypes) = signbit(s1)
@inline function vsignbit(v1::Vec{N,T}) where {N,T<:FloatingTypes}
    U = uint_type(T)
    sm = vbroadcast(Vec{N,U}, sign_mask(T))
    iv = vreinterpret(Vec{N,U}, v1)
    vnot_equal(vand(iv, sm), vbroadcast(Vec{N,U}, 0))
end
@inline function vsignbit(v1::SVec{N,T}) where {N,T<:FloatingTypes}
    SVec(vsignbit(extract_data(v1)))
end
@inline function Base.signbit(v1::SVec{N,T}) where {N,T<:FloatingTypes}
    SVec(vsignbit(extract_data(v1)))
end

@inline vifelse(c::Bool, x, y) = c ? x : y
# @inline vifelse(c::Bool, x, y) = ifelse(c, x, y)
@generated function vifelse(v1::Vec{N,Bool}, v2::Vec{N,T},
        v3::Vec{N,T}) where {N,T}
    btyp = llvmtype(Bool)
    vbtyp = "<$N x $btyp>"
    abtyp = "[$N x $btyp]"
    typ = llvmtype(T)
    vtyp = "<$N x $typ>"
    atyp = "[$N x $typ]"
    decls = String[]
    instrs = String[]
    push!(instrs, "%cond = trunc $vbtyp %0 to <$N x i1>")
    push!(instrs, "%res = select <$N x i1> %cond, $vtyp %1, $vtyp %2")
    push!(instrs, "ret $vtyp %res")
    quote
        $(Expr(:meta, :inline))
        Base.llvmcall($((join(decls, "\n"), join(instrs, "\n"))),
            Vec{N,T},
            Tuple{Vec{N,Bool}, Vec{N,T}, Vec{N,T}},
            v1, v2, v3)
    end
end
@inline function vifelse(v1::AbstractSIMDVector{N,Bool},
        v2::AbstractSIMDVector{N,T}, v3::AbstractSIMDVector{N,T}) where {N,T}
    SVec(vifelse(extract_data(v1), extract_data(v2), extract_data(v3)))
end
# @inline function vifelse(v1::VecOrProd{N,Bool},
#         v2::VecOrProd{N,T}, v3::VecOrProd{N,T}) where {N,T}
#     vifelse(extract_data(v1), extract_data(v2), extract_data(v3))
# end
# @inline function Base.ifelse(v1::SVec{N,Bool},
#         v2::SVec{N,T}, v3::SVec{N,T}) where {N,T}
#     SVec(vifelse(extract_data(v1), extract_data(v2), extract_data(v3)))
# end

@generated function vifelse(mask::U, v2::Vec{N,T}, v3::Vec{N,T}) where {N,T,U<:Unsigned}
    @assert 8sizeof(U) >= N
    btyp = llvmtype(Bool)
    vbtyp = "<$N x $btyp>"
    abtyp = "[$N x $btyp]"
    typ = llvmtype(T)
    vtyp = "<$N x $typ>"
    atyp = "[$N x $typ]"
    decls = String[]
    instrs = String[]

    mtyp_input = llvmtype(U)
    mtyp_trunc = "i$N"

    if mtyp_input == mtyp_trunc
        push!(instrs, "%cond = bitcast $mtyp_input %0 to <$N x i1>")
    else
        push!(instrs, "%condtrunc = trunc $mtyp_input %0 to $mtyp_trunc")
        push!(instrs, "%cond = bitcast $mtyp_trunc %condtrunc to <$N x i1>")
    end
    # push!(instrs, "%cond = trunc $vbtyp %0 to <$N x i1>")
    if T <: FloatingTypes
        push!(instrs, "%res = select fast <$N x i1> %cond, $vtyp %1, $vtyp %2")
    else
        push!(instrs, "%res = select <$N x i1> %cond, $vtyp %1, $vtyp %2")
    end
    push!(instrs, "ret $vtyp %res")
    quote
        $(Expr(:meta, :inline))
        Base.llvmcall($((join(decls, "\n"), join(instrs, "\n"))),
            Vec{$N,$T},
            Tuple{$U, Vec{$N,$T}, Vec{$N,$T}},
            mask, v2, v3)
    end
end
@inline vifelse(U::Unsigned, v2::AbstractSIMDVector, v3::AbstractSIMDVector) = SVec(vifelse(U, extract_data(v2), extract_data(v3)))
@inline vifelse(U::Unsigned, v2::Vec{W,T}, s::Union{T,Int}) where {W,T} = vifelse(U, v2, vbroadcast(Vec{W,T}, s))
@inline vifelse(U::Unsigned, v2::AbstractSIMDVector{W,T}, s::Union{T,Int}) where {W,T} = SVec(vifelse(U, extract_data(v2), vbroadcast(Vec{W,T}, s)))
@inline vifelse(U::Unsigned, s::Union{T,Int}, v2::Vec{W,T}) where {W,T} = vifelse(U, vbroadcast(Vec{W,T}, s), v2)
@inline vifelse(U::Unsigned, s::Union{T,Int}, v2::AbstractSIMDVector{W,T}) where {W,T} = SVec(vifelse(U, vbroadcast(Vec{W,T}, s), extract_data(v2)))

@vectordef visodd function Base.isodd(v) where {N,T<:Integer}
    visequal(vand(v, one(T)), one(T))
end

@inline vand(U::Unsigned, vb::Vec{N,Bool}) where {N} = U & vecbool_to_unsigned(vb)
@inline vand(vb::Vec{N,Bool}, U::Unsigned) where {N} = U & vecbool_to_unsigned(vb)
