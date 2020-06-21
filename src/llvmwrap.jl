# Generic function wrappers

# Functions taking one argument
@inline function llvmwrap(::Val{Op}, v1::_Vec{_W,T1}) where {Op,_W,T1}
    llvmwrap(Val{Op}(), v1, T1)
end
@generated function llvmwrap(::Val{Op}, v1::_Vec{_W,T1}, ::Type{R}) where {Op,_W,T1,R}
    W = _W + 1
    @assert isa(Op, Symbol)
    typ1 = llvmtype(T1)
    vtyp1 = "<$W x $typ1>"
    typr = llvmtype(R)
    vtypr = "<$W x $typr>"
    ins = llvmins(Op, W, T1)
    decls = String[]
    instrs = String[]
    flags = [""]
    if Op in FASTOPS# && T1 <: FloatingTypes
        push!(flags, fastflags(promote_type(T1,R)))
    end
    if ins[1] == '@'
        push!(decls, "declare $vtypr $ins($vtyp1)")
        push!(instrs, "%res = call" * join(flags, " ") * " $vtypr $ins($vtyp1 %0)")
    else
        if Op === :~
            @assert T1 <: IntegerTypes
            otherval = -1
        elseif Op === :inv
            @assert T1 <: FloatingTypes
            otherval = 1.0
        else
            otherval = 0
        end
        otherarg = llvmconst(W, T1, otherval)
        push!(instrs, "%res = $ins" * join(flags, " ") * " $vtyp1 $otherarg, %0")
    end
    push!(instrs, "ret $vtypr %res")
    quote
        $(Expr(:meta, :inline))
        llvmcall(
            $((join(decls, "\n"), join(instrs, "\n"))),
            Vec{$W,$R}, Tuple{Vec{$W,$T1}}, v1
        )
    end
end

# Functions taking one Bool argument
@generated function llvmwrap(::Val{Op}, v1::_Vec{_W,Bool}, ::Type{Bool}) where {Op,_W}
    W = _W + 1
    @assert isa(Op, Symbol)
    btyp = llvmtype(Bool)
    vbtyp = "<$W x $btyp>"
    ins = llvmins(Op, W, Bool)
    decls = String[]
    instrs = String[]
    push!(instrs, "%arg1 = trunc $vbtyp %0 to <$W x i1>")
    otherarg = llvmconst(W, Bool, true)
    push!(instrs, "%res = $ins <$W x i1> $otherarg, %arg1")
    push!(instrs, "%resb = zext <$W x i1> %res to $vbtyp")
    push!(instrs, "ret $vbtyp %resb")
    quote
        $(Expr(:meta, :inline))
        llvmcall($((join(decls, "\n"), join(instrs, "\n"))),
            Vec{$W,Bool}, Tuple{Vec{$W,Bool}}, v1)
    end
end

# Functions taking two arguments
@generated function llvmwrap(
    ::Val{Op}, v1::_Vec{_W,T1}, v2::_Vec{_W,T2}, ::Type{R} = T1
) where {Op,_W,T1,T2,R}
    W = _W + 1
    @assert isa(Op, Symbol)
    typ1 = llvmtype(T1)
    vtyp1 = "<$W x $typ1>"
    typ2 = llvmtype(T2)
    vtyp2 = "<$W x $typ2>"
    typr = llvmtype(R)
    vtypr = "<$W x $typr>"
    ins = llvmins(Op, W, T1)
    decls = String[]
    instrs = String[]
    flags = [""]
    T = promote_type(T1,T2,R)
    if Op in FASTOPS# && T <: FloatingTypes
        push!(flags, fastflags(T))
    end    
    if ins[1] == '@'
        push!(decls, "declare $vtypr $ins($vtyp1, $vtyp2)")
        push!(instrs, "%res = call" * join(flags, " ") * " $vtypr $ins($vtyp1 %0, $vtyp2 %1)")
    else
        push!(instrs, "%res = $ins" * join(flags, " ") * " $vtyp1 %0, %1")
    end
    push!(instrs, "ret $vtypr %res")
    quote
        $(Expr(:meta, :inline))
        llvmcall($((join(decls, "\n"), join(instrs, "\n"))),
            Vec{$W,$R}, Tuple{Vec{$W,$T1}, Vec{$W,$T2}},
            v1, v2)
    end
end
@generated function llvmwrap_notfast(
    ::Val{Op}, v1::_Vec{_W,T1}, v2::_Vec{_W,T2}, ::Type{R} = T1
) where {Op,_W,T1,T2,R}
    @assert isa(Op, Symbol)
    W = _W + 1
    typ1 = llvmtype(T1)
    vtyp1 = "<$W x $typ1>"
    typ2 = llvmtype(T2)
    vtyp2 = "<$W x $typ2>"
    typr = llvmtype(R)
    vtypr = "<$W x $typr>"
    ins = llvmins(Op, W, T1)
    decls = String[]
    instrs = String[]
    if ins[1] == '@'
        push!(decls, "declare $vtypr $ins($vtyp1, $vtyp2)")
        push!(instrs, "%res = call $vtypr $ins($vtyp1 %0, $vtyp2 %1)")
    else
        push!(instrs, "%res = $ins $vtyp1 %0, %1")
    end
    push!(instrs, "ret $vtypr %res")
    quote
        $(Expr(:meta, :inline))
        llvmcall($((join(decls, "\n"), join(instrs, "\n"))),
            Vec{$W,$R}, Tuple{Vec{$W,$T1}, Vec{$W,$T2}},
            v1, v2)
    end
end

# Functions taking two arguments, returning Bool
@generated function llvmwrap(
    ::Val{Op}, v1::_Vec{_W,T1}, v2::_Vec{_W,T2}, ::Type{Bool}
) where {Op,_W,T1,T2}
    @assert isa(Op, Symbol)
    W = _W + 1
    btyp = llvmtype(Bool)
    vbtyp = "<$W x $btyp>"
    typ1 = llvmtype(T1)
    vtyp1 = "<$W x $typ1>"
    typ2 = llvmtype(T2)
    vtyp2 = "<$W x $typ2>"
    ins = llvmins(Op, W, T1)
    decls = String[]
    instrs = String[]
    # if false && W == 1
        # append!(instrs, array2vector("%arg1", W, typ1, "%0", "%arg1arr"))
        # append!(instrs, array2vector("%arg2", W, typ2, "%1", "%arg2arr"))
        # push!(instrs, "%cond = $ins $vtyp1 %arg1, %arg2")
        # push!(instrs, "%res = zext <$W x i1> %cond to $vbtyp")
        # append!(instrs, vector2array("%resarr", W, btyp, "%res"))
        # push!(instrs, "ret $abtyp %resarr")
    # else
    push!(instrs, "%res = $ins $vtyp1 %0, %1")
    push!(instrs, "%resb = zext <$W x i1> %res to $vbtyp")
    push!(instrs, "ret $vbtyp %resb")
    # end
    quote
        $(Expr(:meta, :inline))
        llvmcall($((join(decls, "\n"), join(instrs, "\n"))),
            Vec{$W,Bool}, Tuple{Vec{$W,$T1}, Vec{$W,$T2}},
            v1, v2)
    end
end
# Functions taking two arguments, returning a bitmask
@generated function llvmwrap_bitmask(::Val{Op}, v1::_Vec{_W,T1}, v2::_Vec{_W,T1}) where {Op,_W,T1}
    @assert isa(Op, Symbol)
    W = _W + 1
    typ1 = llvmtype(T1)
    vtyp1 = "<$W x $typ1>"
    ins = llvmins(Op, W, T1)
    decls = String[]
    instrs = String[]
    # if T1 <: FloatingTypes && T2 <: FloatingTypes
        # push!(instrs, "%res = $ins reassoc $vtyp1 %0, %1")
    # else
    push!(instrs, "%res = $ins $vtyp1 %0, %1")
    maskbits = max(8, VectorizationBase.nextpow2(W))
    bitcastname = maskbits == W ? "resu" : "resutrunc"
    push!(instrs, "%$(bitcastname) = bitcast <$W x i1> %res to i$W")
    if maskbits != W
        push!(instrs, "%resu = zext i$W %$(bitcastname) to i$maskbits")
    end
    push!(instrs, "ret i$maskbits %resu")
    julia_mask_type = VectorizationBase.mask_type(maskbits)
    quote
        $(Expr(:meta, :inline))
        llvmcall($((join(decls, "\n"), join(instrs, "\n"))),
            $julia_mask_type, Tuple{Vec{$W,$T1}, Vec{$W,$T1}},
            v1, v2)
    end
end
@generated function llvmwrap_bitmask_notfast(::Val{Op}, v1::_Vec{_W,T1}, v2::_Vec{_W,T1}) where {Op,_W,T1}
    W = _W + 1
    @assert isa(Op, Symbol)
    typ1 = llvmtype(T1)
    vtyp1 = "<$W x $typ1>"
    ins = if Op == :(==)
        "fcmp oeq"
    elseif Op == :(!=)
        "fcmp une"
    elseif Op == :(>)
        "fcmp ogt"
    elseif Op == :(>=)
        "fcmp oge"
    elseif Op == :(<)
        "fcmp olt"
    elseif Op == :(<=)
        "fcmp ole"
    else
        throw("Op $Op not recognized.")
    end
    decls = String[]
    instrs = String[]
    # if T1 <: FloatingTypes && T2 <: FloatingTypes
        # push!(instrs, "%res = $ins reassoc $vtyp1 %0, %1")
    # else
    push!(instrs, "%res = $ins $vtyp1 %0, %1")
    # end
    # push!(instrs, "%resb = zext <$W x i1> %res to $vbtyp")
    maskbits = max(8, VectorizationBase.nextpow2(W))
    bitcastname = maskbits == W ? "resu" : "resutrunc"
    push!(instrs, "%$(bitcastname) = bitcast <$W x i1> %res to i$W")
    if maskbits != W
        push!(instrs, "%resu = zext i$W %$(bitcastname) to i$maskbits")
    end
    push!(instrs, "ret i$maskbits %resu")
    julia_mask_type = VectorizationBase.mask_type(maskbits)
    quote
        $(Expr(:meta, :inline))
        llvmcall($((join(decls, "\n"), join(instrs, "\n"))),
            $julia_mask_type, Tuple{Vec{$W,$T1}, Vec{$W,$T1}},
            v1, v2)
    end
end

# Functions taking two Bool arguments, returning Bool
@generated function llvmwrap(::Val{Op}, v1::_Vec{_W,Bool}, v2::_Vec{_W,Bool}, ::Type{Bool} = Bool) where {Op,_W}
    W = _W + 1
    @assert isa(Op, Symbol)
    btyp = llvmtype(Bool)
    vbtyp = "<$W x $btyp>"
    ins = llvmins(Op, W, Bool)
    decls = String[]
    instrs = String[]
    push!(instrs, "%arg1 = trunc $vbtyp %0 to <$W x i1>")
    push!(instrs, "%arg2 = trunc $vbtyp %1 to <$W x i1>")
    push!(instrs, "%res = $ins <$W x i1> %arg1, %arg2")
    push!(instrs, "%resb = zext <$W x i1> %res to $vbtyp")
    push!(instrs, "ret $vbtyp %resb")
    quote
        $(Expr(:meta, :inline))
        llvmcall($((join(decls, "\n"), join(instrs, "\n"))),
            Vec{$W,Bool}, Tuple{Vec{$W,Bool}, Vec{$W,Bool}},
            v1, v2)
    end
end

# Functions taking three arguments
@generated function llvmwrap(
    ::Val{Op}, v1::_Vec{_W,T1}, v2::_Vec{_W,T2}, v3::_Vec{_W,T3}, ::Type{R} = T1
) where {Op,_W,T1,T2,T3,R}
    W = _W + 1
    @assert isa(Op, Symbol)
    typ1 = llvmtype(T1)
    vtyp1 = "<$W x $typ1>"
    typ2 = llvmtype(T2)
    vtyp2 = "<$W x $typ2>"
    typ3 = llvmtype(T3)
    vtyp3 = "<$W x $typ3>"
    typr = llvmtype(R)
    vtypr = "<$W x $typr>"
    ins = llvmins(Op, W, T1)
    decls = String[]
    instrs = String[]
    if ins[1] == '@'
        push!(decls, "declare $vtypr $ins($vtyp1, $vtyp2, $vtyp3)")
        push!(instrs,
            "%res = call $vtypr $ins($vtyp1 %0, $vtyp2 %1, $vtyp3 %2)")
    else
        push!(instrs, "%res = $ins $vtyp1 %0, %1, %2")
    end
    push!(instrs, "ret $vtypr %res")
    quote
        $(Expr(:meta, :inline))
        llvmcall($((join(decls, "\n"), join(instrs, "\n"))),
            Vec{$W,R},
            Tuple{Vec{$W,$T1}, Vec{$W,$T2}, Vec{$W,$T3}},
            v1, v2, v3)
    end
end
@generated function llvmwrap_fast(
    ::Val{Op}, v1::_Vec{_W,T1}, v2::_Vec{_W,T2}, v3::_Vec{_W,T3}, ::Type{R} = T1
) where {Op,_W,T1,T2,T3,R}
    W = _W + 1
    @assert isa(Op, Symbol)
    typ1 = llvmtype(T1)
    vtyp1 = "<$W x $typ1>"
    typ2 = llvmtype(T2)
    vtyp2 = "<$W x $typ2>"
    typ3 = llvmtype(T3)
    vtyp3 = "<$W x $typ3>"
    typr = llvmtype(R)
    vtypr = "<$W x $typr>"
    ins = llvmins(Op, W, T1)
    decls = String[]
    instrs = String[]
    if ins[1] == '@'
        push!(decls, "declare $vtypr $ins($vtyp1, $vtyp2, $vtyp3)")
        push!(instrs,
            "%res = call fast $vtypr $ins($vtyp1 %0, $vtyp2 %1, $vtyp3 %2)")
    else
        push!(instrs, "%res = $ins fast $vtyp1 %0, %1, %2")
    end
    push!(instrs, "ret $vtypr %res")
    quote
        $(Expr(:meta, :inline))
        llvmcall($((join(decls, "\n"), join(instrs, "\n"))),
            Vec{$W,$R},
            Tuple{Vec{$W,$T1}, Vec{$W,$T2}, Vec{$W,$T3}},
            v1, v2, v3)
    end
end


@generated function llvmwrapshift(
    ::Val{Op}, v1::_Vec{_W,T}, ::Val{I}
) where {Op,_W,T,I}
    W = _W + 1
    @assert isa(Op, Symbol)
    if I >= 0
        op = Op
        i = I
    else
        if Op === :>> || Op === :>>>
            op = :<<
        else
            @assert Op === :<<
            if T <: Unsigned
                op = :>>>
            else
                op = :>>
            end
        end
        i = -I
    end
    @assert op in (:<<, :>>, :>>>)
    @assert i >= 0
    typ = llvmtype(T)
    vtyp = "<$W x $typ>"
    ins = llvmins(op, W, T)
    decls = String[]
    instrs = String[]
    nbits = 8*sizeof(T)
    if (op === :>> && T <: IntTypes) || i < nbits
        count = llvmconst(W, T, min(nbits-1, i))
        push!(instrs, "%res = $ins $vtyp %0, $count")
        push!(instrs, "ret $vtyp %res")
    else
        push!(instrs, "return $vtyp zeroinitializer")
    end
    quote
        $(Expr(:meta, :inline))
        llvmcall($((join(decls, "\n"), join(instrs, "\n"))),
            Vec{$W,$T}, Tuple{Vec{$W,$T}}, v1)
    end
end

@generated function llvmwrapshift(
    ::Val{Op}, v1::_Vec{_W,T}, x2::Unsigned
) where {Op,_W,T}
    W = _W + 1
    @assert isa(Op, Symbol)
    typ = llvmtype(T)
    vtyp = "<$W x $typ>"
    ins = llvmins(Op, W, T)
    decls = String[]
    instrs = String[]
    append!(instrs, scalar2vector("%count", W, typ, "%1"))
    nbits = 8*sizeof(T)
    push!(instrs, "%res = $ins $vtyp %0, %count")
    # push!(instrs, "%inbounds = icmp ult $typ %1, $nbits")
    # if Op === :>> && T <: IntTypes
    #     nbits1 = llvmconst(W, T, 8*sizeof(T)-1)
    #     push!(instrs, "%limit = $ins $vtyp %0, $nbits1")
    #     push!(instrs, "%res = select i1 %inbounds, $vtyp %tmp, $vtyp %limit")
    # else
    #     push!(instrs, "%res = select i1 %inbounds, $vtyp %tmp, $vtyp zeroinitializer")
    # end
    push!(instrs, "ret $vtyp %res")
    quote
        $(Expr(:meta, :inline))
        # Wote that this function might be called with out-of-bounds
        # values for x2, assuming that the results are then ignored
        llvmcall($((join(decls, "\n"), join(instrs, "\n"))),
            Vec{$W,$T}, Tuple{Vec{$W,$T}, $T}, v1, x2 % $T)
    end
end

@generated function llvmwrapshift(
    ::Val{Op}, v1::_Vec{_W,T}, x2::Integer
) where {Op,_W,T}
    W = _W + 1
    if Op === :>> || Op === :>>>
        NegOp = :<<
    else
        @assert Op === :<<
        if T <: Unsigned
            NegOp = :>>>
        else
            NegOp = :>>
        end
    end
    ValOp = Val{Op}()
    ValNegOp = Val{NegOp}()
    quote
        $(Expr(:meta, :inline))
        x2 < 0 ? llvmwrapshift($ValNegOp, v1, unsigned(-x2)) : llvmwrapshift($ValOp, v1, unsigned(x2))
    end
end

@generated function llvmwrapshift(
    ::Val{Op},
    v1::_Vec{_W,T},
    v2::_Vec{_W,U}
) where {Op,_W,T,U<:UIntTypes}
    W = _W + 1
    @assert isa(Op, Symbol)
    typ = llvmtype(T)
    vtyp = "<$W x $typ>"
    ins = llvmins(Op, W, T)
    decls = String[]
    instrs = String[]
    push!(instrs, "%res = $ins $vtyp %0, %1")
    # push!(instrs, "%tmp = $ins $vtyp %0, %1")
    # nbits = llvmconst(W, T, 8*sizeof(T))
    # push!(instrs, "%inbounds = icmp ult $vtyp %1, $nbits")
    # if Op === :>> && T <: IntTypes
    #     nbits1 = llvmconst(W, T, 8*sizeof(T)-1)
    #     push!(instrs, "%limit = $ins $vtyp %0, $nbits1")
    #     push!(instrs,
    #         "%res = select <$W x i1> %inbounds, $vtyp %tmp, $vtyp %limit")
    # else
    #     push!(instrs,
    #         "%res = select <$W x i1> %inbounds, $vtyp %tmp, $vtyp zeroinitializer")
    # end
    push!(instrs, "ret $vtyp %res")
    quote
        $(Expr(:meta, :inline))
        llvmcall($((join(decls, "\n"), join(instrs, "\n"))),
            Vec{$W,$T}, Tuple{Vec{$W,$T}, Vec{$W,$T}},
            v1, vrem(v2, Vec{$W,$T}))
    end
end

@generated function llvmwrapshift(
    ::Val{Op},
    v1::_Vec{_W,T},
    v2::_Vec{_W,U}
) where {Op,_W,T,U<:IntegerTypes}
    W = _W + 1
    if Op === :>> || Op === :>>>
        NegOp = :<<
    else
        @assert Op === :<<
        if T <: Unsigned
            NegOp = :>>>
        else
            NegOp = :>>
        end
    end
    ValOp = Val{Op}()
    ValNegOp = Val{NegOp}()
    quote
        $(Expr(:meta, :inline))
        vifelse(
            vgreater_or_equal(v2, 0),
            llvmwrapshift($ValOp, v1, vrem(v2, Vec{$W,unsigned(U)})),
            llvmwrapshift($ValNegOp, v1, vrem(vsub(v2), Vec{$W,unsigned(U)}))
        )
    end
end


# @generated function llvmwrapshift(
#     ::Val{Op}, v1::_Vec{_W,T}, ::Val{I}
# ) where {Op,W,T,I}
#     @assert isa(Op, Symbol)
#     if I >= 0
#         op = Op
#         i = I
#     else
#         if Op === :>> || Op === :>>>
#             op = :<<
#         else
#             @assert Op === :<<
#             if T <: Unsigned
#                 op = :>>>
#             else
#                 op = :>>
#             end
#         end
#         i = -I
#     end
#     @assert op in (:<<, :>>, :>>>)
#     @assert i >= 0
#     typ = llvmtype(T)
#     vtyp = "<$W x $typ>"
#     ins = llvmins(op, W, T)
#     decls = String[]
#     instrs = String[]
#     nbits = 8*sizeof(T)
#     if (op === :>> && T <: IntTypes) || i < nbits
#         count = llvmconst(W, T, min(nbits-1, i))
#         push!(instrs, "%res = $ins $vtyp %0, $count")
#         push!(instrs, "ret $vtyp %res")
#     else
#         push!(instrs, "return $vtyp zeroinitializer")
#     end
#     quote
#         $(Expr(:meta, :inline))
#         llvmcall($((join(decls, "\n"), join(instrs, "\n"))),
#             Vec{W,T}, Tuple{Vec{W,T}}, v1)
#     end
# end

# @generated function llvmwrapshift(
#     ::Val{Op}, v1::_Vec{_W,T}, x2::Unsigned
# ) where {Op,W,T}
#     @assert isa(Op, Symbol)
#     typ = llvmtype(T)
#     vtyp = "<$W x $typ>"
#     ins = llvmins(Op, W, T)
#     decls = String[]
#     instrs = String[]
#     append!(instrs, scalar2vector("%count", W, typ, "%1"))
#     nbits = 8*sizeof(T)
#     push!(instrs, "%tmp = $ins $vtyp %0, %count")
#     push!(instrs, "%inbounds = icmp ult $typ %1, $nbits")
#     if Op === :>> && T <: IntTypes
#         nbits1 = llvmconst(W, T, 8*sizeof(T)-1)
#         push!(instrs, "%limit = $ins $vtyp %0, $nbits1")
#         push!(instrs, "%res = select i1 %inbounds, $vtyp %tmp, $vtyp %limit")
#     else
#         push!(instrs, "%res = select i1 %inbounds, $vtyp %tmp, $vtyp zeroinitializer")
#     end
#     push!(instrs, "ret $vtyp %res")
#     quote
#         $(Expr(:meta, :inline))
#         # Note that this function might be called with out-of-bounds
#         # values for x2, assuming that the results are then ignored
#         llvmcall($((join(decls, "\n"), join(instrs, "\n"))),
#             Vec{W,T}, Tuple{Vec{W,T}, T}, v1, x2 % T)
#     end
# end

# @generated function llvmwrapshift(
#     ::Val{Op}, v1::_Vec{_W,T}, x2::Integer
# ) where {Op,W,T}
#     if Op === :>> || Op === :>>>
#         NegOp = :<<
#     else
#         @assert Op === :<<
#         if T <: Unsigned
#             NegOp = :>>>
#         else
#             NegOp = :>>
#         end
#     end
#     ValOp = Val{Op}()
#     ValNegOp = Val{NegOp}()
#     quote
#         $(Expr(:meta, :inline))
#         x2 < 0 ? llvmwrapshift($ValNegOp, v1, unsigned(-x2)) : llvmwrapshift($ValOp, v1, unsigned(x2))
#     end
# end

# @generated function llvmwrapshift(
#     ::Val{Op},
#     v1::_Vec{_W,T},
#     v2::_Vec{_W,U}
# ) where {Op,W,T,U<:UIntTypes}
#     @assert isa(Op, Symbol)
#     typ = llvmtype(T)
#     vtyp = "<$W x $typ>"
#     ins = llvmins(Op, W, T)
#     decls = String[]
#     instrs = String[]
#     push!(instrs, "%res = $ins $vtyp %0, %1")
#     # push!(instrs, "%tmp = $ins $vtyp %0, %1")
#     # nbits = llvmconst(W, T, 8*sizeof(T))
#     # push!(instrs, "%inbounds = icmp ult $vtyp %1, $nbits")
#     # if Op === :>> && T <: IntTypes
#     #     nbits1 = llvmconst(W, T, 8*sizeof(T)-1)
#     #     push!(instrs, "%limit = $ins $vtyp %0, $nbits1")
#     #     push!(instrs,
#     #         "%res = select <$W x i1> %inbounds, $vtyp %tmp, $vtyp %limit")
#     # else
#     #     push!(instrs,
#     #         "%res = select <$W x i1> %inbounds, $vtyp %tmp, $vtyp zeroinitializer")
#     # end
#     push!(instrs, "ret $vtyp %res")
#     quote
#         $(Expr(:meta, :inline))
#         llvmcall($((join(decls, "\n"), join(instrs, "\n"))),
#             Vec{W,T}, Tuple{Vec{W,T}, Vec{W,T}},
#             v1, vrem(v2, Vec{W,T}))
#     end
# end

# @generated function llvmwrapshift(
#     ::Val{Op},
#     v1::_Vec{_W,T},
#     v2::_Vec{_W,U}
# ) where {Op,W,T,U<:IntegerTypes}
#     if Op === :>> || Op === :>>>
#         NegOp = :<<
#     else
#         @assert Op === :<<
#         if T <: Unsigned
#             NegOp = :>>>
#         else
#             NegOp = :>>
#         end
#     end
#     ValOp = Val{Op}()
#     ValNegOp = Val{NegOp}()
#     quote
#         $(Expr(:meta, :inline))
#         vifelse(
#             vgreater_or_equal(v2, 0),
#             llvmwrapshift($ValOp, v1, vrem(v2, Vec{W,unsigned(U)})),
#             llvmwrapshift($ValNegOp, v1, vrem(vsub(v2), Vec{W,unsigned(U)}))
#         )
#     end
# end


