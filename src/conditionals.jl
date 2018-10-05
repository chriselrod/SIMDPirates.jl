
# Conditionals

for op ∈ (:(==), :(!=), :(<), :(<=), :(>), :(>=) )
    rename = VECTOR_SYMBOLS[op]
    @eval begin
        # scalar versions handled in floating_point_arithmetic.jl
        # @inline $rename(s1::ScalarTypes, s2::ScalarTypes) = $op(s1,s2)
        @inline $rename(v1::Vec{N,T}, v2::Vec{N,T}) where {N,T} =
            llvmwrap(Val{$(QuoteNode(op))}, v1, v2, Bool)
    end
end
@inline visfinite(s::ScalarTypes) = isfinite(s)
@inline function visfinite(v1::Vec{N,T}) where {N,T<:FloatingTypes}
    U = uint_type(T)
    em = Vec{N,U}(exponent_mask(T))
    iv = pirate_reinterpret(Vec{N,U}, v1)
    vnot_equal(vand(iv, em), em)
end
@inline visinf(s1::ScalarTypes) = isinf(s1)
@inline visinf(v1::Vec{N,T}) where {N,T<:FloatingTypes} = vequal(vabs(v1), vbroadcast(Vec{N,T},Inf))
@inline visnan(s1::ScalarTypes) = isnan(s1)
@inline visnan(v1::Vec{N,T}) where {N,T<:FloatingTypes} = vnot_equal(v1, v1)
@inline vissubnormal(s1::ScalarTypes) = issubnormal(s1)
@inline function vissubnormal(v1::Vec{N,T}) where {N,T<:FloatingTypes}
    U = uint_type(T)
    em = Vec{N,U}(exponent_mask(T))
    sm = Vec{N,U}(significand_mask(T))
    iv = pirate_reinterpret(Vec{N,U}, v1)
    vand(vequal(vand(iv, em), vbroadcast(Vec{N,U}, 0)), vnot_equal(vand(iv, sm), vbroadcast(Vec{N,U}, 0)))
end
@inline signbit(s1::ScalarTypes) = signbit(s1)
@inline function vsignbit(v1::Vec{N,T}) where {N,T<:FloatingTypes}
    U = uint_type(T)
    sm = Vec{N,U}(sign_mask(T))
    iv = pirate_reinterpret(Vec{N,U}, v1)
    vnot_equal(vand(iv, sm), vbroadcast(Vec{N,U}, 0))
end


@inline vifelse(c::Bool, x, y) = ifelse(c, x, y)
@generated function vifelse(v1::Vec{N,Bool}, v2::Vec{N,T},
        v3::Vec{N,T}) where {N,T}
    btyp = llvmtype(Bool)
    vbtyp = "<$N x $btyp>"
    abtyp = "[$N x $btyp]"
    typ = llvmtype(T)
    vtyp = "<$N x $typ>"
    atyp = "[$N x $typ]"
    decls = []
    instrs = []
    if false && N == 1
        append!(instrs, array2vector("%arg1", N, btyp, "%0", "%arg1arr"))
        append!(instrs, array2vector("%arg2", N, typ, "%1", "%arg2arr"))
        append!(instrs, array2vector("%arg3", N, typ, "%2", "%arg3arr"))
        push!(instrs, "%cond = trunc $vbtyp %arg1 to <$N x i1>")
        push!(instrs, "%res = select <$N x i1> %cond, $vtyp %arg2, $vtyp %arg3")
        append!(instrs, vector2array("%resarr", N, typ, "%res"))
        push!(instrs, "ret $atyp %resarr")
    else
        push!(instrs, "%cond = trunc $vbtyp %0 to <$N x i1>")
        push!(instrs, "%res = select <$N x i1> %cond, $vtyp %1, $vtyp %2")
        push!(instrs, "ret $vtyp %res")
    end
    quote
        $(Expr(:meta, :inline))
        Base.llvmcall($((join(decls, "\n"), join(instrs, "\n"))),
            Vec{N,T},
            Tuple{Vec{N,Bool}, Vec{N,T}, Vec{N,T}},
            v1, v2, v3)
    end
end