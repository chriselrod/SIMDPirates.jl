# Generic function wrappers

# Functions taking one argument
@generated function llvmwrap(::Type{Val{Op}}, v1::Vec{N,T1},
        ::Type{R} = T1) where {Op,N,T1,R}
    @assert isa(Op, Symbol)
    typ1 = llvmtype(T1)
    vtyp1 = "<$N x $typ1>"
    typr = llvmtype(R)
    vtypr = "<$N x $typr>"
    ins = llvmins(Val{Op}, N, T1)
    decls = []
    instrs = []
    if ins[1] == '@'
        push!(decls, "declare $vtypr $ins($vtyp1)")
        push!(instrs, "%res = call $vtypr $ins($vtyp1 %0)")
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
        otherarg = llvmconst(N, T1, otherval)
        push!(instrs, "%res = $ins $vtyp1 $otherarg, %0")
    end
    push!(instrs, "ret $vtypr %res")
    quote
        $(Expr(:meta, :inline))
        Base.llvmcall($((join(decls, "\n"), join(instrs, "\n"))),
            Vec{N,R}, Tuple{Vec{N,T1}}, v1)
    end
end

# Functions taking one Bool argument
@generated function llvmwrap(::Type{Val{Op}}, v1::Vec{N,Bool},
        ::Type{Bool} = Bool) where {Op,N}
    @assert isa(Op, Symbol)
    btyp = llvmtype(Bool)
    vbtyp = "<$N x $btyp>"
    ins = llvmins(Val{Op}, N, Bool)
    decls = []
    instrs = []
    push!(instrs, "%arg1 = trunc $vbtyp %0 to <$N x i1>")
    otherarg = llvmconst(N, Bool, true)
    push!(instrs, "%res = $ins <$N x i1> $otherarg, %arg1")
    push!(instrs, "%resb = zext <$N x i1> %res to $vbtyp")
    push!(instrs, "ret $vbtyp %resb")
    quote
        $(Expr(:meta, :inline))
        Base.llvmcall($((join(decls, "\n"), join(instrs, "\n"))),
            Vec{N,Bool}, Tuple{Vec{N,Bool}}, v1)
    end
end

# Functions taking two arguments
@generated function llvmwrap(::Type{Val{Op}}, v1::Vec{N,T1},
        v2::Vec{N,T2}, ::Type{R} = T1) where {Op,N,T1,T2,R}
    @assert isa(Op, Symbol)
    typ1 = llvmtype(T1)
    vtyp1 = "<$N x $typ1>"
    typ2 = llvmtype(T2)
    vtyp2 = "<$N x $typ2>"
    typr = llvmtype(R)
    vtypr = "<$N x $typr>"
    ins = llvmins(Val{Op}, N, T1)
    decls = []
    instrs = []
    if ins[1] == '@'
        push!(decls, "declare $vtypr $ins($vtyp1, $vtyp2)")
        push!(instrs, "%res = call $vtypr $ins($vtyp1 %0, $vtyp2 %1)")
    else
        push!(instrs, "%res = $ins $vtyp1 %0, %1")
    end
    push!(instrs, "ret $vtypr %res")
    quote
        $(Expr(:meta, :inline))
        Base.llvmcall($((join(decls, "\n"), join(instrs, "\n"))),
            Vec{N,R}, Tuple{Vec{N,T1}, Vec{N,T2}},
            v1, v2)
    end
end

# Functions taking two arguments, returning Bool
@generated function llvmwrap(::Type{Val{Op}}, v1::Vec{N,T1},
        v2::Vec{N,T2}, ::Type{Bool}) where {Op,N,T1,T2}
    @assert isa(Op, Symbol)
    btyp = llvmtype(Bool)
    vbtyp = "<$N x $btyp>"
    abtyp = "[$N x $btyp]"
    typ1 = llvmtype(T1)
    vtyp1 = "<$N x $typ1>"
    atyp1 = "[$N x $typ1]"
    typ2 = llvmtype(T2)
    vtyp2 = "<$N x $typ2>"
    atyp2 = "[$N x $typ2]"
    ins = llvmins(Val{Op}, N, T1)
    decls = []
    instrs = []
    if false && N == 1
        append!(instrs, array2vector("%arg1", N, typ1, "%0", "%arg1arr"))
        append!(instrs, array2vector("%arg2", N, typ2, "%1", "%arg2arr"))
        push!(instrs, "%cond = $ins $vtyp1 %arg1, %arg2")
        push!(instrs, "%res = zext <$N x i1> %cond to $vbtyp")
        append!(instrs, vector2array("%resarr", N, btyp, "%res"))
        push!(instrs, "ret $abtyp %resarr")
    else
        push!(instrs, "%res = $ins $vtyp1 %0, %1")
        push!(instrs, "%resb = zext <$N x i1> %res to $vbtyp")
        push!(instrs, "ret $vbtyp %resb")
    end
    quote
        $(Expr(:meta, :inline))
        Base.llvmcall($((join(decls, "\n"), join(instrs, "\n"))),
            Vec{N,Bool}, Tuple{Vec{N,T1}, Vec{N,T2}},
            v1, v2)
    end
end

# Functions taking two Bool arguments, returning Bool
@generated function llvmwrap(::Type{Val{Op}}, v1::Vec{N,Bool},
        v2::Vec{N,Bool}, ::Type{Bool} = Bool) where {Op,N}
    @assert isa(Op, Symbol)
    btyp = llvmtype(Bool)
    vbtyp = "<$N x $btyp>"
    ins = llvmins(Val{Op}, N, Bool)
    decls = []
    instrs = []
    push!(instrs, "%arg1 = trunc $vbtyp %0 to <$N x i1>")
    push!(instrs, "%arg2 = trunc $vbtyp %1 to <$N x i1>")
    push!(instrs, "%res = $ins <$N x i1> %arg1, %arg2")
    push!(instrs, "%resb = zext <$N x i1> %res to $vbtyp")
    push!(instrs, "ret $vbtyp %resb")
    quote
        $(Expr(:meta, :inline))
        Base.llvmcall($((join(decls, "\n"), join(instrs, "\n"))),
            Vec{N,Bool}, Tuple{Vec{N,Bool}, Vec{N,Bool}},
            v1, v2)
    end
end

# Functions taking three arguments
@generated function llvmwrap(::Type{Val{Op}}, v1::Vec{N,T1},
        v2::Vec{N,T2}, v3::Vec{N,T3}, ::Type{R} = T1) where {Op,N,T1,T2,T3,R}
    @assert isa(Op, Symbol)
    typ1 = llvmtype(T1)
    vtyp1 = "<$N x $typ1>"
    typ2 = llvmtype(T2)
    vtyp2 = "<$N x $typ2>"
    typ3 = llvmtype(T3)
    vtyp3 = "<$N x $typ3>"
    typr = llvmtype(R)
    vtypr = "<$N x $typr>"
    ins = llvmins(Val{Op}, N, T1)
    decls = []
    instrs = []
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
        Base.llvmcall($((join(decls, "\n"), join(instrs, "\n"))),
            Vec{N,R},
            Tuple{Vec{N,T1}, Vec{N,T2}, Vec{N,T3}},
            v1, v2, v3)
    end
end

@generated function llvmwrapshift(::Type{Val{Op}}, v1::Vec{N,T},
                                  ::Type{Val{I}}) where {Op,N,T,I}
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
    vtyp = "<$N x $typ>"
    ins = llvmins(Val{op}, N, T)
    decls = []
    instrs = []
    nbits = 8*sizeof(T)
    if (op === :>> && T <: IntTypes) || i < nbits
        count = llvmconst(N, T, min(nbits-1, i))
        push!(instrs, "%res = $ins $vtyp %0, $count")
        push!(instrs, "ret $vtyp %res")
    else
        zero = llvmconst(N, T, 0)
        push!(instrs, "return $vtyp $zero")
    end
    quote
        $(Expr(:meta, :inline))
        Base.llvmcall($((join(decls, "\n"), join(instrs, "\n"))),
            Vec{N,T}, Tuple{Vec{N,T}}, v1)
    end
end

@generated function llvmwrapshift(::Type{Val{Op}}, v1::Vec{N,T},
                                  x2::Unsigned) where {Op,N,T}
    @assert isa(Op, Symbol)
    typ = llvmtype(T)
    vtyp = "<$N x $typ>"
    ins = llvmins(Val{Op}, N, T)
    decls = []
    instrs = []
    append!(instrs, scalar2vector("%count", N, typ, "%1"))
    nbits = 8*sizeof(T)
    push!(instrs, "%tmp = $ins $vtyp %0, %count")
    push!(instrs, "%inbounds = icmp ult $typ %1, $nbits")
    if Op === :>> && T <: IntTypes
        nbits1 = llvmconst(N, T, 8*sizeof(T)-1)
        push!(instrs, "%limit = $ins $vtyp %0, $nbits1")
        push!(instrs, "%res = select i1 %inbounds, $vtyp %tmp, $vtyp %limit")
    else
        zero = llvmconst(N, T, 0)
        push!(instrs, "%res = select i1 %inbounds, $vtyp %tmp, $vtyp $zero")
    end
    push!(instrs, "ret $vtyp %res")
    quote
        $(Expr(:meta, :inline))
        # Note that this function might be called with out-of-bounds
        # values for x2, assuming that the results are then ignored
        Base.llvmcall($((join(decls, "\n"), join(instrs, "\n"))),
            Vec{N,T}, Tuple{Vec{N,T}, T}, v1, x2 % T)
    end
end

@generated function llvmwrapshift(::Type{Val{Op}}, v1::Vec{N,T},
                                  x2::Integer) where {Op,N,T}
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
    ValOp = Val{Op}
    ValNegOp = Val{NegOp}
    quote
        $(Expr(:meta, :inline))
        ifelse(x2 >= 0,
               llvmwrapshift($ValOp, v1, unsigned(x2)),
               llvmwrapshift($ValNegOp, v1, unsigned(-x2)))
    end
end

@generated function llvmwrapshift(::Type{Val{Op}},
                                  v1::Vec{N,T},
                                  v2::Vec{N,U}) where {Op,N,T,U<:UIntTypes}
    @assert isa(Op, Symbol)
    typ = llvmtype(T)
    vtyp = "<$N x $typ>"
    ins = llvmins(Val{Op}, N, T)
    decls = []
    instrs = []
    push!(instrs, "%tmp = $ins $vtyp %0, %1")
    nbits = llvmconst(N, T, 8*sizeof(T))
    push!(instrs, "%inbounds = icmp ult $vtyp %1, $nbits")
    if Op === :>> && T <: IntTypes
        nbits1 = llvmconst(N, T, 8*sizeof(T)-1)
        push!(instrs, "%limit = $ins $vtyp %0, $nbits1")
        push!(instrs,
            "%res = select <$N x i1> %inbounds, $vtyp %tmp, $vtyp %limit")
    else
        zero = llvmconst(N, T, 0)
        push!(instrs,
            "%res = select <$N x i1> %inbounds, $vtyp %tmp, $vtyp $zero")
    end
    push!(instrs, "ret $vtyp %res")
    quote
        $(Expr(:meta, :inline))
        Base.llvmcall($((join(decls, "\n"), join(instrs, "\n"))),
            Vec{N,T}, Tuple{Vec{N,T}, Vec{N,T}},
            v1, vrem(v2, Vec{N,T}))
    end
end

@generated function llvmwrapshift(::Type{Val{Op}},
                                  v1::Vec{N,T},
                                  v2::Vec{N,U}) where {Op,N,T,U<:IntegerTypes}
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
    ValOp = Val{Op}
    ValNegOp = Val{NegOp}
    quote
        $(Expr(:meta, :inline))
        vifelse(vgreater_or_equal(v2, 0),
                llvmwrapshift($ValOp, v1, v2 % Vec{N,unsigned(U)}),
                llvmwrapshift($ValNegOp, v1, -v2 % Vec{N,unsigned(U)}))
    end
end
