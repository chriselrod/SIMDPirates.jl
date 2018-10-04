
# Conditionals

for op in (:(==), :(!=), :(<), :(<=), :(>), :(>=))
    @eval begin
        @inline $op(v1::Vec{N,T}, v2::Vec{N,T}) where {N,T} =
            llvmwrap(Val{$(QuoteNode(op))}, v1, v2, Bool)
    end
end
@inline function isfinite(v1::Vec{N,T}) where {N,T<:FloatingTypes}
    U = uint_type(T)
    em = Vec{N,U}(exponent_mask(T))
    iv = pirate_reinterpret(Vec{N,U}, v1)
    iv & em != em
end
@inline isinf(v1::Vec{N,T}) where {N,T<:FloatingTypes} = abs(v1) == broadcast(Vec{N,T},Inf)
@inline isnan(v1::Vec{N,T}) where {N,T<:FloatingTypes} = v1 != v1
@inline function issubnormal(v1::Vec{N,T}) where {N,T<:FloatingTypes}
    U = uint_type(T)
    em = Vec{N,U}(exponent_mask(T))
    sm = Vec{N,U}(significand_mask(T))
    iv = pirate_reinterpret(Vec{N,U}, v1)
    (iv & em == broadcast(Vec{N,U}, 0)) & (iv & sm != broadcast(Vec{N,U}, 0))
end
@inline function signbit(v1::Vec{N,T}) where {N,T<:FloatingTypes}
    U = uint_type(T)
    sm = Vec{N,U}(sign_mask(T))
    iv = pirate_reinterpret(Vec{N,U}, v1)
    iv & sm != broadcast(Vec{N,U}, 0)
end


vifelse(c::Bool, x, y) = ifelse(c, x, y)
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
