


"""
I'm not sure on the details, but I think this function can only allocate up to

524288 doubles

that is it may not be able to allocate arrays with more than about half a million elements, or 4 megabytes.

Another limitation is that LLVM will assume that each `alloca` within a function can be merged.
Therefore,
ptr1 = alloca(100)
ptr2 = alloca(2_000)
ptr3 = alloca(1_000)
ptr4 = alloca(500)

will result in allocating 2000 doubles, and ptr1 == ptr2 == ptr3 == ptr4.
This is very likely not what someone writing the above intended, but I do not yet know a workaround.
"""
@generated function alloca(::Val{N}, ::Type{T} = Float64, ::Val{Align} = Val{64}()) where {N, T, Align}
    typ = llvmtype(T)
    ptyp = JuliaPointerType
    instrs = String[]
    push!(instrs, "%ptr = alloca $typ, i32 $N, align $Align")
    push!(instrs, "%iptr = ptrtoint $typ* %ptr to $ptyp")
    # push!(instrs, "%ptr = alloca i8, i32 $(N*sizeof(T)), align $Align")
    # push!(instrs, "%iptr = ptrtoint i8* %ptr to $ptyp")
    push!(instrs, "ret $ptyp %iptr")
    quote
        $(Expr(:meta,:inline))
        Base.llvmcall(
            $(join(instrs, "\n")),
            Ptr{$T}, Tuple{}
        )
    end
end
@generated function alloca(N::Int32, ::Type{T} = Float64, ::Val{Align} = Val{64}()) where {T, Align}
    typ = llvmtype(T)
    ptyp = JuliaPointerType
    instrs = String[]
    push!(instrs, "%ptr = alloca $typ, i32 %0, align $Align")
    push!(instrs, "%iptr = ptrtoint $typ* %ptr to $ptyp")
    # push!(instrs, "%ptr = alloca i8, i32 $(N*sizeof(T)), align $Align")
    # push!(instrs, "%iptr = ptrtoint i8* %ptr to $ptyp")
    push!(instrs, "ret $ptyp %iptr")
    quote
        $(Expr(:meta,:inline))
        Base.llvmcall(
            $(join(instrs, "\n")),
            Ptr{$T}, Tuple{Int32}, N
        )
    end
end
@inline function alloca(N::Integer, ::Type{T} = Float64, ::Val{Align} = Val{64}()) where {T, Align}
    alloca(N % Int32, T, Val{Align}())
end

"""
ReadOrWrite: 0 for read, 1 for write
"""
@generated function prefetch(ptr::Ptr{T}, ::Val{Locality} = Val(1), ::Val{ReadOrWrite} = Val(0)) where {T, Locality, ReadOrWrite}
    prefetch_call_string = """%addr = inttoptr i$(8sizeof(Int)) %0 to i8*
    call void @llvm.prefetch(i8* %addr, i32 $ReadOrWrite, i32 $Locality, i32 1)
    ret void"""
    quote
        $(Expr(:meta, :inline))
        Base.llvmcall(
            ("declare void @llvm.prefetch(i8*, i32, i32, i32)",
             $prefetch_call_string), Cvoid, Tuple{Ptr{$T}}, ptr
        )
    end
end
@inline function prefetch(ptr::Union{VectorizationBase.AbstractStridedPointer,Ptr}, i, ::Val{Locality} = Val(1), ::Val{ReadOrWrite} = Val(0)) where {Locality, ReadOrWrite}
    prefetch(gep(ptr, i), Val{Locality}(), Val{ReadOrWrite}())
end


function valloc(::Type{T}, N::Int) where {T}
    @assert N > 0
    W, Wshift = VectorizationBase.pick_vector_width_shift(T)
    # We use padding to align the address of the first element, and
    # also to ensure that we can access past the last element up to
    # the next full vector width
    mem = Vector{T}(undef, (N + W) & (-W))
    addr = reinterpret(Int, pointer(mem))
    reg_size = VectorizationBase.REGISTER_SIZE
    off = (addr & (reg_size - 1)) >>> VectorizationBase.intlog2(sizeof(T))
    off = W - off
    view(mem, off + 1 : off + N)
end

# @generated function sub2ind(dims::NTuple{N}, I::NTuple{N}) where {N}
#     ex = :(I[$N] - 1)
#     for i = (N - 1):-1:1
#         ex = :(I[$i] - 1 + dims[$i] * $ex)
#     end
#     quote
#         $(Expr(:meta, :inline))
#         $ex + 1
#     end
# end
@generated function vload(
    ::Type{Vec{W,T}}, ptr::Ptr{T}, ::Val{Aligned}, ::Val{Nontemporal}
) where {W,T,Aligned, Nontemporal}
    @assert isa(Aligned, Bool)
    ptyp = JuliaPointerType
    typ = llvmtype(T)
    vtyp = "<$W x $typ>"
    decl = """
    !1 = !{!\"noaliasdomain\"}
    !2 = !{!\"noaliasscope\", !1}
    !3 = !{!2}
    """
    # decls = String[]
    instrs = String[]
    if Aligned
        align = Base.datatype_alignment(Vec{W,T})
    else
        align = sizeof(T)   # This is overly optimistic
    end
    flags = [""]
    align > 0 && push!(flags, "align $align")
    Nontemporal && push!(flags, "!nontemporal !{i32 1}")
    push!(flags, "!alias.scope !3")
    push!(instrs, "%ptr = inttoptr $ptyp %0 to $vtyp*")
    push!(instrs, "%res = load $vtyp, $vtyp* %ptr" * join(flags, ", "))
    push!(instrs, "ret $vtyp %res")
    quote
        $(Expr(:meta, :inline))
        Base.llvmcall($((decl, join(instrs, "\n"))),
            Vec{$W,$T}, Tuple{Ptr{$T}}, ptr)
    end
end
@generated function vload(
    ptr::Ptr{T}, i::_MM{W}, ::Val{Aligned}, ::Val{Nontemporal}
) where {W,T,Aligned, Nontemporal}
    @assert isa(Aligned, Bool)
    ptyp = JuliaPointerType
    typ = llvmtype(T)
    vtyp = "<$W x $typ>"
    decl = """
    !1 = !{!\"noaliasdomain\"}
    !2 = !{!\"noaliasscope\", !1}
    !3 = !{!2}
    """
    instrs = String[]
    if Aligned
        align = Base.datatype_alignment(Vec{W,T})
    else
        align = sizeof(T)   # This is overly optimistic
    end
    flags = [""]
    align > 0 && push!(flags, "align $align")
    Nontemporal && push!(flags, "!nontemporal !{i32 1}")
    push!(flags, "!alias.scope !3")
    push!(instrs, "%typptr = inttoptr $ptyp %0 to $typ*")
    push!(instrs, "%offsetptr = getelementptr inbounds $typ, $typ* %typptr, $ptyp %1")
    push!(instrs, "%ptr = bitcast $typ* %offsetptr to $vtyp*")
    # push!(instrs, "%offsetptr = add nsw nuw $ptyp %0, %1") # uncommenting requires *8n
    # push!(instrs, "%ptr = inttoptr $ptyp %offsetptr to $vtyp*")
    push!(instrs, "%res = load $vtyp, $vtyp* %ptr" * join(flags, ", "))
    push!(instrs, "ret $vtyp %res")
    quote
        $(Expr(:meta, :inline))
        SVec(Base.llvmcall($((decl, join(instrs, "\n"))),
            Vec{$W,$T}, Tuple{Ptr{$T},Int}, ptr, (i.i % Int)))
    end
end
@generated function vload(
    ::Type{Vec{W,T}}, ptr::Ptr{T}, mask::U, ::Val{Aligned}
) where {W,T,Aligned,U<:Unsigned}
    @assert isa(Aligned, Bool)
    @assert 8sizeof(U) >= W
    ptyp = JuliaPointerType
    typ = llvmtype(T)
    vtyp = "<$W x $typ>"
    mtyp_input = llvmtype(U)
    mtyp_trunc = "i$W"
    decls = String["!1 = !{!\"noaliasdomain\"}","!2 = !{!\"noaliasscope\", !1}", "!3 = !{!2}"]
    # decls = String[]
    instrs = String[]
    if Aligned
        align = Base.datatype_alignment(Vec{W,T})
    else
        align = sizeof(T)   # This is overly optimistic
    end
    push!(instrs, "%ptr = inttoptr $ptyp %0 to $vtyp*")
    if mtyp_input == mtyp_trunc
        push!(instrs, "%mask = bitcast $mtyp_input %1 to <$W x i1>")
    else
        push!(instrs, "%masktrunc = trunc $mtyp_input %1 to $mtyp_trunc")
        push!(instrs, "%mask = bitcast $mtyp_trunc %masktrunc to <$W x i1>")
    end
    push!(decls,
        "declare $vtyp @llvm.masked.load.$(suffix(W,T))($vtyp*, i32, <$W x i1>, $vtyp)"
    )
    push!(instrs,
        "%res = call $vtyp @llvm.masked.load.$(suffix(W,T))($vtyp* %ptr, i32 $align, <$W x i1> %mask, $vtyp zeroinitializer), !alias.scope !3"#undef)"# 
    )
    push!(instrs, "ret $vtyp %res")
    quote
        $(Expr(:meta, :inline))
        Base.llvmcall($((join(decls, "\n"), join(instrs, "\n"))),
            Vec{$W,$T}, Tuple{Ptr{$T}, $U}, ptr, mask)
    end
end
@generated function vload(
    ptr::Ptr{T}, i::_MM{W}, mask::U, ::Val{Aligned}
) where {W,T,U<:Unsigned,Aligned}
    ptyp = JuliaPointerType
    typ = llvmtype(T)
    vtyp = "<$W x $typ>"
    mtyp_input = llvmtype(U)
    mtyp_trunc = "i$W"
    decls = String["!1 = !{!\"noaliasdomain\"}","!2 = !{!\"noaliasscope\", !1}", "!3 = !{!2}"]
    # decls = String[]
    instrs = String[]
    if Aligned
        align = Base.datatype_alignment(Vec{W,T})
    else
        align = sizeof(T)   # This is overly optimistic
    end
    push!(instrs, "%typptr = inttoptr $ptyp %0 to $typ*")
    push!(instrs, "%offsetptr = getelementptr inbounds $typ, $typ* %typptr, $ptyp %1")
    push!(instrs, "%ptr = bitcast $typ* %offsetptr to $vtyp*")
    if mtyp_input == mtyp_trunc
        push!(instrs, "%mask = bitcast $mtyp_input %2 to <$W x i1>")
    else
        push!(instrs, "%masktrunc = trunc $mtyp_input %2 to $mtyp_trunc")
        push!(instrs, "%mask = bitcast $mtyp_trunc %masktrunc to <$W x i1>")
    end
    push!(decls, "declare $vtyp @llvm.masked.load.$(suffix(W,T))($vtyp*, i32, <$W x i1>, $vtyp)")
    push!(instrs,"%res = call $vtyp @llvm.masked.load.$(suffix(W,T))($vtyp* %ptr, i32 $align, <$W x i1> %mask, $vtyp zeroinitializer), !alias.scope !3")
    push!(instrs, "ret $vtyp %res")
    quote
        $(Expr(:meta, :inline))
        SVec(Base.llvmcall($((join(decls, "\n"), join(instrs, "\n"))),
            Vec{$W,$T}, Tuple{Ptr{$T}, Int, $U}, ptr, i.i, mask))
    end
end

for v ∈ (:Vec, :SVec, :Val, :_MM)
    pargs = if v === :Val
        Union{Symbol,Expr}[ :(::Val{W}), :(ptr::Ptr{T})]
    elseif v === :_MM
        Union{Symbol,Expr}[:(ptr::Ptr{T})]
    else
        Union{Symbol,Expr}[:(::Type{$v{W,T}}), :(ptr::Ptr{T})]
    end
    for index ∈ (true,false)
        if v === :_MM
            index || continue
            icall = Union{Symbol,Expr}[:ptr, :i]#Expr(:call, :gep, :ptr, :i)]
            iargs = push!(copy(pargs), :(i::Union{_MM{W},AbstractSIMDVector{W,<:Integer}}))
        else
            if index
                icall = Union{Symbol,Expr}[:ptr, :(_MM{W}(i))]#Expr(:call, :gep, :ptr, :i)]
                iargs = push!(copy(pargs), :i)
            else
                icall = Union{Symbol,Expr}[:(Vec{W,T}), :ptr]
                iargs = pargs
            end
        end
        for mask ∈ [ :nothing, :Unsigned, :Mask ]
            if mask !== :nothing
                margs = push!(copy(iargs), Expr(:(::), :mask, mask))
                if mask === :Mask
                    mcall = push!(copy(icall), Expr(:call, :extract_data, :mask))
                else
                    mcall = push!(copy(icall), :mask)
                end
                fopts = [:vload,:vloada]
            else
                margs = iargs
                mcall = icall
                fopts = [:vload,:vloada,:vloadnt]
            end
            for f ∈ fopts
                body = Expr(:call, :vload, mcall...)
                push!(body.args, Expr(:call, Expr(:curly, :Val, f !== :vload))) # aligned arg
                if mask === :nothing
                    push!(body.args, Expr(:call, Expr(:curly, :Val, f === :vloadnt))) # nontemporal argf
                end
                if index
                    if v === :Vec
                        body = Expr(:call, :extract_data, body)
                    end
                else
                    if v !== :Vec
                        body = Expr(:call, :SVec, body)
                    end
                end
                @eval @inline function $f($(margs...)) where {W,T}
                    $body
                end
            end
        end
    end
end

@generated function vstore!(
    ptr::Ptr{T}, v::Vec{W,T}, ::Val{Aligned}, ::Val{Nontemporal}
) where {W,T,Aligned, Nontemporal}
    @assert isa(Aligned, Bool)
    ptyp = JuliaPointerType
    typ = llvmtype(T)
    vtyp = "<$W x $typ>"
    instrs = String[]
    if Aligned# || Nontemporal
        align = Base.datatype_alignment(Vec{W,T})
    else
        align = sizeof(T)   # This is overly optimistic
    end
    flags = [""]
    align > 0 && push!(flags, "align $align")
    Nontemporal && push!(flags, "!nontemporal !{i32 1}")
    push!(instrs, "%ptr = inttoptr $ptyp %0 to $vtyp*")
    push!(instrs, "store $vtyp %1, $vtyp* %ptr" * join(flags, ", "))
    push!(instrs, "ret void")
    quote
        $(Expr(:meta, :inline))
        Base.llvmcall(
            $(join(instrs, "\n")),
            Cvoid, Tuple{Ptr{$T}, Vec{$W,$T}}, ptr, v
        )
    end
end
@generated function vstore!(
    ptr::Ptr{T}, v::Vec{W,T}, i::I, ::Val{Aligned}, ::Val{Nontemporal}
) where {W,T,I,Aligned, Nontemporal}
    @assert isa(Aligned, Bool)
    ityp = llvmtype(I)
    ptyp = JuliaPointerType
    typ = llvmtype(T)
    vtyp = "<$W x $typ>"
    instrs = String[]
    if Aligned# || Nontemporal
        align = Base.datatype_alignment(Vec{W,T})
    else
        align = sizeof(T)   # This is overly optimistic
    end
    flags = [""]
    align > 0 && push!(flags, "align $align")
    Nontemporal && push!(flags, "!nontemporal !{i32 1}")
    push!(instrs, "%typptr = inttoptr $ptyp %0 to $typ*")
    push!(instrs, "%offsetptr = getelementptr inbounds $typ, $typ* %typptr, $ityp %2")
    push!(instrs, "%ptr = bitcast $typ* %offsetptr to $vtyp*")
    push!(instrs, "store $vtyp %1, $vtyp* %ptr" * join(flags, ", "))
    push!(instrs, "ret void")
    quote
        $(Expr(:meta, :inline))
        Base.llvmcall(
            $(join(instrs, "\n")),
            Cvoid, Tuple{Ptr{$T}, Vec{$W,$T}, $I}, ptr, v, i
        )
    end
end

@generated function vstore!(
    ptr::Ptr{T}, v::Vec{W,T}, mask::U, ::Val{Aligned}# = Val{false}()
) where {W,T,Aligned,U<:Unsigned}
    @assert isa(Aligned, Bool)
    ptyp = JuliaPointerType
    typ = llvmtype(T)
    vtyp = "<$W x $typ>"
    mtyp_input = llvmtype(U)
    mtyp_trunc = "i$W"
    instrs = String[]
    if Aligned
        align = Base.datatype_alignment(Vec{W,T})
    else
        align = sizeof(T)   # This is overly optimistic
    end
    push!(instrs, "%ptr = inttoptr $ptyp %0 to $vtyp*")
    if mtyp_input == mtyp_trunc
        push!(instrs, "%mask = bitcast $mtyp_input %2 to <$W x i1>")
    else
        push!(instrs, "%masktrunc = trunc $mtyp_input %2 to $mtyp_trunc")
        push!(instrs, "%mask = bitcast $mtyp_trunc %masktrunc to <$W x i1>")
    end
    decl = "declare void @llvm.masked.store.$(suffix(W,T))($vtyp, $vtyp*, i32, <$W x i1>)"
    push!(instrs,
        "call void @llvm.masked.store.$(suffix(W,T))($vtyp %1, $vtyp* %ptr, i32 $align, <$W x i1> %mask)"
    )
    push!(instrs, "ret void")
    quote
        $(Expr(:meta, :inline))
        Base.llvmcall($((decl, join(instrs, "\n"))),
            Cvoid, Tuple{Ptr{$T}, Vec{$W,$T}, $U},
            ptr, v, mask)
    end
end
@generated function vstore!(
    ptr::Ptr{T}, v::Vec{W,T}, i::I, mask::U, ::Val{Aligned}# = Val{false}()
) where {W,T,Aligned,U<:Unsigned,I<:Integer}
    @assert isa(Aligned, Bool)
    ityp = llvmtype(I)
    ptyp = JuliaPointerType
    typ = llvmtype(T)
    vtyp = "<$W x $typ>"
    mtyp_input = llvmtype(U)
    mtyp_trunc = "i$W"
    instrs = String[]
    if Aligned
        align = Base.datatype_alignment(Vec{W,T})
    else
        align = sizeof(T)   # This is overly optimistic
    end
    push!(instrs, "%typptr = inttoptr $ptyp %0 to $typ*")
    push!(instrs, "%offsetptr = getelementptr inbounds $typ, $typ* %typptr, $ityp %2")
    push!(instrs, "%ptr = bitcast $typ* %offsetptr to $vtyp*")
    if mtyp_input == mtyp_trunc
        push!(instrs, "%mask = bitcast $mtyp_input %3 to <$W x i1>")
    else
        push!(instrs, "%masktrunc = trunc $mtyp_input %3 to $mtyp_trunc")
        push!(instrs, "%mask = bitcast $mtyp_trunc %masktrunc to <$W x i1>")
    end
    decl = "declare void @llvm.masked.store.$(suffix(W,T))($vtyp, $vtyp*, i32, <$W x i1>)"
    push!(instrs,
        "call void @llvm.masked.store.$(suffix(W,T))($vtyp %1, $vtyp* %ptr, i32 $align, <$W x i1> %mask)"
    )
    push!(instrs, "ret void")
    quote
        $(Expr(:meta, :inline))
        Base.llvmcall($((decl, join(instrs, "\n"))),
            Cvoid, Tuple{Ptr{$T}, Vec{$W,$T}, $I, $U},
            ptr, v, i, mask)
    end
end

@generated function vnoaliasstore!(
    ptr::Ptr{T}, v::Vec{W,T}, i::I, ::Val{Aligned}, ::Val{Nontemporal}
) where {W,T,I,Aligned, Nontemporal}
    @assert isa(Aligned, Bool)
    ityp = llvmtype(I)
    ptyp = JuliaPointerType
    typ = llvmtype(T)
    vtyp = "<$W x $typ>"
    decls = String["!1 = !{!\"noaliasdomain\"}","!2 = !{!\"noaliasscope\", !1}", "!3 = !{!2}"]
    instrs = String[]
    if Aligned# || Nontemporal
        align = Base.datatype_alignment(Vec{W,T})
    else
        align = sizeof(T)   # This is overly optimistic
    end
    flags = [""]
    align > 0 && push!(flags, "align $align")
    Nontemporal && push!(flags, "!nontemporal !{i32 1}")
    push!(flags, "!noalias !3")
    push!(instrs, "%typptr = inttoptr $ptyp %0 to $typ*")
    push!(instrs, "%offsetptr = getelementptr inbounds $typ, $typ* %typptr, $ityp %2")
    push!(instrs, "%ptr = bitcast $typ* %offsetptr to $vtyp*")
    push!(instrs, "store $vtyp %1, $vtyp* %ptr" * join(flags, ", "))
    push!(instrs, "ret void")
    quote
        $(Expr(:meta, :inline))
        Base.llvmcall(
            $((join(decls, "\n"), join(instrs, "\n"))),
            Cvoid, Tuple{Ptr{$T}, Vec{$W,$T}, $I}, ptr, v, i
        )
    end
end

@generated function vnoaliasstore!(
    ptr::Ptr{T}, v::Vec{W,T}, i::I, mask::U, ::Val{Aligned}# = Val{false}()
) where {W,T,Aligned,U<:Unsigned,I<:Integer}
    @assert isa(Aligned, Bool)
    ityp = llvmtype(I)
    ptyp = JuliaPointerType
    typ = llvmtype(T)
    vtyp = "<$W x $typ>"
    mtyp_input = llvmtype(U)
    mtyp_trunc = "i$W"
    decls = String["!1 = !{!\"noaliasdomain\"}","!2 = !{!\"noaliasscope\", !1}", "!3 = !{!2}"]
    instrs = String[]
    if Aligned
        align = Base.datatype_alignment(Vec{W,T})
    else
        align = sizeof(T)   # This is overly optimistic
    end
    # push!(flags, "!noalias !3")
    push!(instrs, "%typptr = inttoptr $ptyp %0 to $typ*")
    push!(instrs, "%offsetptr = getelementptr inbounds $typ, $typ* %typptr, $ityp %2")
    push!(instrs, "%ptr = bitcast $typ* %offsetptr to $vtyp*")
    if mtyp_input == mtyp_trunc
        push!(instrs, "%mask = bitcast $mtyp_input %3 to <$W x i1>")
    else
        push!(instrs, "%masktrunc = trunc $mtyp_input %3 to $mtyp_trunc")
        push!(instrs, "%mask = bitcast $mtyp_trunc %masktrunc to <$W x i1>")
    end
    push!(decls,
        "declare void @llvm.masked.store.$(suffix(W,T))($vtyp, $vtyp*, i32, <$W x i1>)"
    )
    push!(instrs,
        "call void @llvm.masked.store.$(suffix(W,T))($vtyp %1, $vtyp* %ptr, i32 $align, <$W x i1> %mask), !noalias !3"
    )
    push!(instrs, "ret void")
    quote
        $(Expr(:meta, :inline))
        Base.llvmcall($((join(decls, "\n"), join(instrs, "\n"))),
            Cvoid, Tuple{Ptr{$T}, Vec{$W,$T}, $I, $U},
            ptr, v, i, mask)
    end
end


let pargs = Union{Symbol,Expr}[:(ptr::Ptr{T}), :(v::AbstractSIMDVector{W,T})]
    for index ∈ 0:2
        if index == 0
            icall = Union{Symbol,Expr}[:ptr, :(extract_data(v))]
            iargs = pargs
        elseif index == 1
            icall = Union{Symbol,Expr}[:ptr, :(extract_data(v)), :i]
            iargs = push!(copy(pargs), :(i::Integer))
        else#if index == 2
            icall = Union{Symbol,Expr}[:ptr, :(extract_data(v)), :(i.i)]
            iargs = push!(copy(pargs), :(i::_MM{W}))
        # else
            # icall = Union{Symbol,Expr}[:ptr, :(extract_data(v)), Expr(:(.), :i, QuoteNode(:i))]
            # iargs = push!(copy(pargs), Expr(:(::), :i, Expr(:curly, :_MM, :W)))
        end
        for mask ∈ [ :nothing, :Unsigned, :Mask ]
            for prefix ∈ (:vstore, :vnoaliasstore)
                if mask !== :nothing
                    margs = push!(copy(iargs), Expr(:(::), :mask, mask))
                    if mask === :Mask
                        mcall = push!(copy(icall), Expr(:call, :extract_data, :mask))
                    else
                        mcall = push!(copy(icall), :mask)
                    end
                    suffix = [:!,:a!]
                else
                    margs = iargs
                    mcall = icall
                    suffix = [:!,:a!,:nt!]
                end
                for s ∈ suffix
                    f = Symbol(prefix, s)
                    body = Expr(:call, Symbol(prefix, :!))
                    append!(body.args, mcall)
                    push!(body.args, Expr(:call, Expr(:curly, :Val, s !== :!))) # aligned arg
                    if mask === :nothing
                        push!(body.args, Expr(:call, Expr(:curly, :Val, s === :nt!))) # nontemporal argf
                    end
                    @eval @inline function $f($(margs...)) where {W,T}
                        $body
                    end
                end
            end
        end
    end
end

# @generated function vloadscope(
#     ::Type{Vec{W,T}}, ptr::Ptr{T}, ::Val{Scope},
#     ::Val{Aligned}, ::Val{Nontemporal}
# ) where {W, T, Scope, Aligned, Nontemporal}
#     @assert isa(Aligned, Bool)
#     ptyp = JuliaPointerType
#     typ = llvmtype(T)
#     vtyp = "<$W x $typ>"
#     domain,scope,list = Scope::NTuple{3,Int}
#     decls = String["!$domain = !{!$domain}","!$scope = !{!$scope, !$domain}", "!$list = !{!$scope}"]
#     # decls = String[]
#     instrs = String[]
#     if Aligned
#         align = Base.datatype_alignment(Vec{W,T})
#     else
#         align = sizeof(T)   # This is overly optimistic
#     end
#     flags = [""]
#     align > 0 && push!(flags, "align $align")
#     push!(flags, "!alias.scope !$list")
#     Nontemporal && push!(flags, "!nontemporal !{i32 1}")
#     push!(instrs, "%ptr = inttoptr $ptyp %0 to $vtyp*")
#     push!(instrs, "%res = load $vtyp, $vtyp* %ptr" * join(flags, ", "))
#     push!(instrs, "ret $vtyp %res")
#     quote
#         $(Expr(:meta, :inline))
#         Base.llvmcall($((join(decls, "\n"), join(instrs, "\n"))),
#             Vec{$W,$T}, Tuple{Ptr{$T}}, ptr)
#     end
# end


@generated function gather(
   ptr::Vec{W,Ptr{T}}, ::Val{Aligned}# = Val{false}()
) where {W,T,Aligned}
    @assert isa(Aligned, Bool)
    ptyp = JuliaPointerType
    vptyp = "<$W x $ptyp>"
    typ = llvmtype(T)
    vtyp = "<$W x $typ>"
    vptrtyp = "<$W x $typ*>"
    instrs = String[]
    if Aligned
        align = Base.datatype_alignment(Vec{W,T})
    else
        align = sizeof(T)   # This is overly optimistic
    end
    mask = join((", i1 true" for i ∈ 2:W))
    push!(instrs, "%ptr = inttoptr $vptyp %0 to $vptrtyp")
    decl = "declare $vtyp @llvm.masked.gather.$(suffix(W,T)).$(suffix(W,Ptr{T}))($vptrtyp, i32, <$W x i1>, $vtyp)"
    push!(instrs, "%res = call $vtyp @llvm.masked.gather.$(suffix(W,T)).$(suffix(W,Ptr{T}))($vptrtyp %ptr, i32 $align, <$W x i1> <i1 true$(mask)>, $vtyp undef)")
    push!(instrs, "ret $vtyp %res")
    quote
        $(Expr(:meta, :inline))
        Base.llvmcall(
            $((decl, join(instrs, "\n"))),
            Vec{$W,$T}, Tuple{Vec{$W,Ptr{$T}}}, ptr
        )
    end
end

@generated function gather(
   ptr::Vec{W,Ptr{T}}, mask::U, ::Val{Aligned}# = Val{false}()
) where {W,T,Aligned,U<:Unsigned}
    @assert isa(Aligned, Bool)
    @assert 8sizeof(U) >= W
    ptyp = JuliaPointerType
    vptyp = "<$W x $ptyp>"
    typ = llvmtype(T)
    vtyp = "<$W x $typ>"
    vptrtyp = "<$W x $typ*>"
    mtyp_input = llvmtype(U)
    mtyp_trunc = "i$W"
    instrs = String[]
    if Aligned
        align = Base.datatype_alignment(Vec{W,T})
    else
        align = sizeof(T)   # This is overly optimistic
    end
    push!(instrs, "%ptr = inttoptr $vptyp %0 to $vptrtyp")
    if mtyp_input == mtyp_trunc
        push!(instrs, "%mask = bitcast $mtyp_input %1 to <$W x i1>")
    else
        push!(instrs, "%masktrunc = trunc $mtyp_input %1 to $mtyp_trunc")
        push!(instrs, "%mask = bitcast $mtyp_trunc %masktrunc to <$W x i1>")
    end
    decls = "declare $vtyp @llvm.masked.gather.$(suffix(W,T)).$(suffix(W,Ptr{T}))($vptrtyp, i32, <$W x i1>, $vtyp)"
    push!(instrs, "%res = call $vtyp @llvm.masked.gather.$(suffix(W,T)).$(suffix(W,Ptr{T}))($vptrtyp %ptr, i32 $align, <$W x i1> %mask, $vtyp zeroinitializer)")#undef)")
    push!(instrs, "ret $vtyp %res")
    quote
        $(Expr(:meta, :inline))
        Base.llvmcall(
            $((decl, join(instrs, "\n"))),
            Vec{$W,$T}, Tuple{Vec{$W,Ptr{$T}}, $U}, ptr, mask
        )
    end
end
@generated function vload(
   ptr::Ptr{T}, i::Vec{W,I}, ::Val{Aligned}# = Val{false}()
) where {W,T,I<:Integer,Aligned}
    @assert isa(Aligned, Bool)
    ptyp = JuliaPointerType
    vptyp = "<$W x $ptyp>"
    typ = llvmtype(T)
    vtyp = "<$W x $typ>"
    vptrtyp = "<$W x $typ*>"
    instrs = String[]
    if Aligned
        align = Base.datatype_alignment(Vec{W,T})
    else
        align = sizeof(T)   # This is overly optimistic
    end
    ityp = llvmtype(I)
    vityp = "<$W x $ityp>"
    push!(instrs, "%sptr = inttoptr $ptyp %0 to $typ*")
    push!(instrs, "%ptr = getelementptr inbounds $typ, $typ* %sptr, $vityp %1")
    mask = join((", i1 true" for i ∈ 2:W))
    decl = "declare $vtyp @llvm.masked.gather.$(suffix(W,T)).$(suffix(W,Ptr{T}))($vptrtyp, i32, <$W x i1>, $vtyp)"
    push!(instrs, "%res = call $vtyp @llvm.masked.gather.$(suffix(W,T)).$(suffix(W,Ptr{T}))($vptrtyp %ptr, i32 $align, <$W x i1> <i1 true$(mask)>, $vtyp undef)")
    push!(instrs, "ret $vtyp %res")
    quote
        $(Expr(:meta, :inline))
        Base.llvmcall(
            $((decl, join(instrs, "\n"))),
            Vec{$W,$T}, Tuple{Ptr{$T},Vec{$W,$I}}, ptr, i
        )
    end
end
@inline vload(ptr::Ptr, i::SVec{W,I}, ::Val{Aligned}) where {W,I<:Integer,Aligned} = SVec(vload(ptr, extract_data(i), Val{Aligned}()))
@inline vload(ptr::Ptr, i::SVec{W,I}, mask::Unsigned, ::Val{Aligned}) where {W,I<:Integer,Aligned} = SVec(vload(ptr, extract_data(i), mask, Val{Aligned}()))
@inline vload(ptr::Ptr, i::SVec{W,I}, mask::Mask{W}, ::Val{Aligned}) where {W,I<:Integer,Aligned} = SVec(vload(ptr, extract_data(i), mask.u, Val{Aligned}()))
@inline vload(ptr::Ptr, i::SVec{W,I}, ::Val{Aligned}, ::Val{false}) where {W,I<:Integer,Aligned} = SVec(vload(ptr, extract_data(i), Val{Aligned}()))
@inline vload(ptr::Ptr, i::SVec{W,I}, mask::Unsigned, ::Val{Aligned}, ::Val{false}) where {W,I<:Integer,Aligned} = SVec(vload(ptr, extract_data(i), mask, Val{Aligned}()))
@inline vload(ptr::Ptr, i::SVec{W,I}, mask::Mask{W}, ::Val{Aligned}, ::Val{false}) where {W,I<:Integer,Aligned} = SVec(vload(ptr, extract_data(i), mask.u, Val{Aligned}()))
@generated function vload(
   ptr::Ptr{T}, i::Vec{W,I}, mask::U, ::Val{Aligned}# = Val{false}()
) where {W,T,I<:Integer,Aligned,U<:Unsigned}
    @assert isa(Aligned, Bool)
    @assert 8sizeof(U) >= W
    ptyp = JuliaPointerType
    vptyp = "<$W x $ptyp>"
    typ = llvmtype(T)
    vtyp = "<$W x $typ>"
    vptrtyp = "<$W x $typ*>"
    mtyp_input = llvmtype(U)
    mtyp_trunc = "i$W"
    instrs = String[]
    if Aligned
        align = Base.datatype_alignment(Vec{W,T})
    else
        align = sizeof(T)   # This is overly optimistic
    end
    ityp = llvmtype(I)
    vityp = "<$W x $ityp>"
    push!(instrs, "%sptr = inttoptr $ptyp %0 to $typ*")
    push!(instrs, "%ptr = getelementptr inbounds $typ, $typ* %sptr, $vityp %1")
    if mtyp_input == mtyp_trunc
        push!(instrs, "%mask = bitcast $mtyp_input %2 to <$W x i1>")
    else
        push!(instrs, "%masktrunc = trunc $mtyp_input %2 to $mtyp_trunc")
        push!(instrs, "%mask = bitcast $mtyp_trunc %masktrunc to <$W x i1>")
    end
    decl = "declare $vtyp @llvm.masked.gather.$(suffix(W,T)).$(suffix(W,Ptr{T}))($vptrtyp, i32, <$W x i1>, $vtyp)"
    push!(instrs, "%res = call $vtyp @llvm.masked.gather.$(suffix(W,T)).$(suffix(W,Ptr{T}))($vptrtyp %ptr, i32 $align, <$W x i1> %mask, $vtyp zeroinitializer)")#undef)")
    push!(instrs, "ret $vtyp %res")
    quote
        $(Expr(:meta, :inline))
        Base.llvmcall(
            $((decl, join(instrs, "\n"))),
            Vec{$W,$T}, Tuple{Ptr{$T},Vec{$W,$I},$U}, ptr, i, mask
        )
    end
end
@generated function scatter!(
    ptr::Vec{W,Ptr{T}}, v::Vec{W,T}, ::Val{Aligned}# = Val{false}()
) where {W,T,Aligned}
    @assert isa(Aligned, Bool)
    ptyp = JuliaPointerType
    vptyp = "<$W x $ptyp>"
    typ = llvmtype(T)
    vtyp = "<$W x $typ>"
    vptrtyp = "<$W x $typ*>"
    instrs = String[]
    if Aligned
        align = Base.datatype_alignment(Vec{W,T})
    else
        align = sizeof(T)   # This is overly optimistic
    end
    push!(instrs, "%ptr = inttoptr $vptyp %1 to $vptrtyp")
    mask = join((", i1 true" for i ∈ 2:W))
    # push!(decls, "declare void @llvm.masked.scatter.$(suffix(W,T)).$(suffix(W,Ptr{T}))($vtyp, $vptrtyp, i32, <$W x i1>)")
    # push!(instrs, "call void @llvm.masked.scatter.$(suffix(W,T)).$(suffix(W,Ptr{T}))($vtyp %0, $vptrtyp %ptr, i32 $align, <$W x i1> <i1 true$(mask)>)")
    decl = "declare void @llvm.masked.scatter.$(suffix(W,T))($vtyp, $vptrtyp, i32, <$W x i1>)"
    push!(instrs, "call void @llvm.masked.scatter.$(suffix(W,T))($vtyp %0, $vptrtyp %ptr, i32 $align, <$W x i1> <i1 true$(mask)>)")
    push!(instrs, "ret void")
    quote
        $(Expr(:meta, :inline))
        Base.llvmcall($((decl, join(instrs, "\n"))),
            Cvoid, Tuple{Vec{$W,$T}, Vec{$W,Ptr{$T}}}, v, ptr)
    end
end
@generated function scatter!(
    ptr::Vec{W,Ptr{T}}, v::Vec{W,T}, mask::U, ::Val{Aligned}# = Val{false}()
) where {W,T,Aligned,U<:Unsigned}
    @assert isa(Aligned, Bool)
    @assert 8sizeof(U) >= W
    ptyp = JuliaPointerType
    vptyp = "<$W x $ptyp>"
    typ = llvmtype(T)
    vtyp = "<$W x $typ>"
    vptrtyp = "<$W x $typ*>"
    mtyp_input = llvmtype(U)
    mtyp_trunc = "i$W"
    instrs = String[]
    if Aligned
        align = Base.datatype_alignment(Vec{W,T})
    else
        align = sizeof(T)   # This is overly optimistic
    end
    push!(instrs, "%ptr = inttoptr $vptyp %1 to $vptrtyp")
    if mtyp_input == mtyp_trunc
        push!(instrs, "%mask = bitcast $mtyp_input %2 to <$W x i1>")
    else
        push!(instrs, "%masktrunc = trunc $mtyp_input %2 to $mtyp_trunc")
        push!(instrs, "%mask = bitcast $mtyp_trunc %masktrunc to <$W x i1>")
    end
    # push!(decls,
        # "declare void @llvm.masked.scatter.$(suffix(W,T)).$(suffix(W,Ptr{T}))($vtyp, $vptrtyp, i32, <$W x i1>)")
    # push!(instrs,
        # "call void @llvm.masked.scatter.$(suffix(W,T)).$(suffix(W,Ptr{T}))($vtyp %0, $vptrtyp %ptr, i32 $align, <$W x i1> %mask)")
    decl = "declare void @llvm.masked.scatter.$(suffix(W,T))($vtyp, $vptrtyp, i32, <$W x i1>)"
    push!(instrs, "call void @llvm.masked.scatter.$(suffix(W,T))($vtyp %0, $vptrtyp %ptr, i32 $align, <$W x i1> %mask)")
    push!(instrs, "ret void")
    quote
        $(Expr(:meta, :inline))
        Base.llvmcall(
            $((decl, join(instrs, "\n"))),
            Cvoid, Tuple{Vec{$W,$T}, Vec{$W,Ptr{$T}}, $U}, v, ptr, mask
        )
    end
end
@generated function vstore!(
    ptr::Ptr{T}, v::Vec{W,T}, i::Vec{W,I}, ::Val{Aligned}# = Val{false}()
) where {W,T,Aligned, I<:Integer}
    @assert isa(Aligned, Bool)
    ptyp = JuliaPointerType
    vptyp = "<$W x $ptyp>"
    typ = llvmtype(T)
    vtyp = "<$W x $typ>"
    vptrtyp = "<$W x $typ*>"
    instrs = String[]
    if Aligned
        align = Base.datatype_alignment(Vec{W,T})
    else
        align = sizeof(T)   # This is overly optimistic
    end
    ityp = llvmtype(I)
    vityp = "<$W x $ityp>"
    push!(instrs, "%sptr = inttoptr $ptyp %0 to $typ*")
    push!(instrs, "%ptr = getelementptr inbounds $typ, $typ* %sptr, $vityp %2")
    mask = join((", i1 true" for i ∈ 2:W))
    # push!(decls, "declare void @llvm.masked.scatter.$(suffix(W,T)).$(suffix(W,Ptr{T}))($vtyp, $vptrtyp, i32, <$W x i1>)")
    # push!(instrs, "call void @llvm.masked.scatter.$(suffix(W,T)).$(suffix(W,Ptr{T}))($vtyp %0, $vptrtyp %ptr, i32 $align, <$W x i1> <i1 true$(mask)>)")
    decl = "declare void @llvm.masked.scatter.$(suffix(W,T))($vtyp, $vptrtyp, i32, <$W x i1>)"
    push!(instrs, "call void @llvm.masked.scatter.$(suffix(W,T))($vtyp %1, $vptrtyp %ptr, i32 $align, <$W x i1> <i1 true$(mask)>)")
    push!(instrs, "ret void")
    quote
        $(Expr(:meta, :inline))
        Base.llvmcall($((decl, join(instrs, "\n"))),
            Cvoid, Tuple{Ptr{$T}, Vec{$W,$T}, Vec{$W,$I}}, ptr, v, i)
    end
end
@generated function vstore!(
    ptr::Ptr{T}, v::Vec{W,T}, i::Vec{W,I}, mask::U, ::Val{Aligned}# = Val{false}()
) where {W,T,Aligned,U<:Unsigned,I<:Integer}
    @assert isa(Aligned, Bool)
    @assert 8sizeof(U) >= W
    ptyp = JuliaPointerType
    vptyp = "<$W x $ptyp>"
    typ = llvmtype(T)
    vtyp = "<$W x $typ>"
    vptrtyp = "<$W x $typ*>"
    mtyp_input = llvmtype(U)
    mtyp_trunc = "i$W"
    instrs = String[]
    if Aligned
        align = Base.datatype_alignment(Vec{W,T})
    else
        align = sizeof(T)   # This is overly optimistic
    end
    ityp = llvmtype(I)
    vityp = "<$W x $ityp>"
    push!(instrs, "%sptr = inttoptr $ptyp %0 to $typ*")
    push!(instrs, "%ptr = getelementptr inbounds $typ, $typ* %sptr, $vityp %2")
    if mtyp_input == mtyp_trunc
        push!(instrs, "%mask = bitcast $mtyp_input %3 to <$W x i1>")
    else
        push!(instrs, "%masktrunc = trunc $mtyp_input %3 to $mtyp_trunc")
        push!(instrs, "%mask = bitcast $mtyp_trunc %masktrunc to <$W x i1>")
    end
    # push!(decls,
        # "declare void @llvm.masked.scatter.$(suffix(W,T)).$(suffix(W,Ptr{T}))($vtyp, $vptrtyp, i32, <$W x i1>)")
    # push!(instrs,
        # "call void @llvm.masked.scatter.$(suffix(W,T)).$(suffix(W,Ptr{T}))($vtyp %0, $vptrtyp %ptr, i32 $align, <$W x i1> %mask)")
    decl = "declare void @llvm.masked.scatter.$(suffix(W,T))($vtyp, $vptrtyp, i32, <$W x i1>)"
    push!(instrs, "call void @llvm.masked.scatter.$(suffix(W,T))($vtyp %1, $vptrtyp %ptr, i32 $align, <$W x i1> %mask)")
    push!(instrs, "ret void")
    quote
        $(Expr(:meta, :inline))
        Base.llvmcall(
            $((decl, join(instrs, "\n"))),
            Cvoid, Tuple{Ptr{$T}, Vec{$W,$T}, Vec{$W,$I}, $U}, ptr, v, i, mask
        )
    end
end
@inline vstore!(ptr::Ptr{T}, v::AbstractSIMDVector{W,T}, i::AbstractSIMDVector{W,I}) where {W,T,I<:Integer} = vstore!(ptr, extract_data(v), extract_data(i), Val{false}())
@inline vstore!(ptr::Ptr{T}, v::AbstractSIMDVector{W,T}, i::AbstractSIMDVector{W,I}, mask::Mask{W}) where {W,T,I<:Integer} = vstore!(ptr, extract_data(v), extract_data(i), mask.u, Val{false}())
@inline vstore!(ptr::Ptr{T}, v::AbstractSIMDVector{W,T}, i::AbstractSIMDVector{W,I}, mask::Unsigned) where {W,T,I<:Integer} = vstore!(ptr, extract_data(v), extract_data(i), mask, Val{false}())

@inline vstore!(ptr::Ptr{T}, v::AbstractSIMDVector{W,T}, i::AbstractSIMDVector{W,I}, ::Val{Aligned}) where {W,T,I<:Integer,Aligned} = vstore!(ptr, extract_data(v), extract_data(i), Val{Aligned}())
@inline vstore!(ptr::Ptr{T}, v::AbstractSIMDVector{W,T}, i::AbstractSIMDVector{W,I}, mask::Mask{W}, ::Val{Aligned}) where {W,T,I<:Integer,Aligned} = vstore!(ptr, extract_data(v), extract_data(i), mask.u, Val{Aligned}())
@inline vstore!(ptr::Ptr{T}, v::AbstractSIMDVector{W,T}, i::AbstractSIMDVector{W,I}, mask::Unsigned, ::Val{Aligned}) where {W,T,I<:Integer,Aligned} = vstore!(ptr, extract_data(v), extract_data(i), mask, Val{Aligned}())

@inline vstore!(ptr::Ptr{T}, v::AbstractSIMDVector{W,T}, i::AbstractSIMDVector{W,I}, ::Val{Aligned}, ::Val{false}) where {W,T,I<:Integer,Aligned} = vstore!(ptr, extract_data(v), extract_data(i), Val{Aligned}())
@inline vstore!(ptr::Ptr{T}, v::AbstractSIMDVector{W,T}, i::AbstractSIMDVector{W,I}, mask::Mask{W}, ::Val{Aligned}, ::Val{false}) where {W,T,I<:Integer,Aligned} = vstore!(ptr, extract_data(v), extract_data(i), mask.u, Val{Aligned}())
@inline vstore!(ptr::Ptr{T}, v::AbstractSIMDVector{W,T}, i::AbstractSIMDVector{W,I}, mask::Unsigned, ::Val{Aligned}, ::Val{false}) where {W,T,I<:Integer,Aligned} = vstore!(ptr, extract_data(v), extract_data(i), mask, Val{Aligned}())

@generated function lifetime_start!(ptr::Ptr{T}, ::Val{L}) where {L,T}
    ptyp = JuliaPointerType
    decl = "declare void @llvm.lifetime.start(i64, i8* nocapture)"
    instrs = [
        "%ptr = inttoptr $ptyp %0 to i8*",
        "call void @llvm.lifetime.start(i64 $(L*sizeof(T)), i8* %ptr)",
        "ret void"
    ]
    quote
        $(Expr(:meta,:inline))
        Base.llvmcall(
            $((decl, join(instrs, "\n"))),
            Cvoid, Tuple{Ptr{$T}}, ptr
        )
    end
end
@generated function lifetime_end!(ptr::Ptr{T}, ::Val{L}) where {L,T}
    ptyp = JuliaPointerType
    decl = "declare void @llvm.lifetime.end(i64, i8* nocapture)"
    instrs = [
        "%ptr = inttoptr $ptyp %0 to i8*",
        "call void @llvm.lifetime.end(i64 $(L*sizeof(T)), i8* %ptr)",
        "ret void"
    ]
    quote
        $(Expr(:meta,:inline))
        Base.llvmcall(
            $((decl, join(instrs, "\n"))),
            Cvoid, Tuple{Ptr{$T}}, ptr
        )
    end
end

@inline function lifetime_start!(ptr::AbstractPointer{T}) where {T}
    lifetime_start!(pointer(ptr), Val{1}())
end
@inline function lifetime_start!(ptr::AbstractPointer{T}, ::Val{L}) where {T,L}
    lifetime_start!(pointer(ptr), Val{L}())
end
@inline function lifetime_start!(ptr::Ptr{T}) where {T}
    lifetime_start!(ptr, Val{1}())
end

@inline function lifetime_end!(ptr::AbstractPointer{T}) where {T}
   lifetime_end!(pointer(ptr), Val{1}())
end
@inline function lifetime_end!(ptr::AbstractPointer{T}, ::Val{L}) where {T,L}
    lifetime_end!(pointer(ptr), Val{L}())
end
@inline function lifetime_end!(ptr::Ptr{T}) where {T}
    lifetime_end!(ptr, Val{1}())
end
# Fallback is to do nothing
@inline lifetime_start!(::Any) = nothing
@inline lifetime_end!(::Any) = nothing

@generated function compressstore!(
    ptr::Ptr{T}, v::Vec{W,T}, mask::U
) where {W,T,U<:Unsigned}
    @assert 8sizeof(U) >= W
    ptyp = JuliaPointerType
    typ = llvmtype(T)
    vtyp = "<$W x $typ>"
    mtyp_input = llvmtype(U)
    mtyp_trunc = "i$W"
    instrs = String[]
    push!(instrs, "%ptr = inttoptr $ptyp %1 to $typ*")
    if mtyp_input == mtyp_trunc
        push!(instrs, "%mask = bitcast $mtyp_input %2 to <$W x i1>")
    else
        push!(instrs, "%masktrunc = trunc $mtyp_input %2 to $mtyp_trunc")
        push!(instrs, "%mask = bitcast $mtyp_trunc %masktrunc to <$W x i1>")
    end
    decl = "declare void @llvm.masked.compressstore.$(suffix(W,T))($vtyp, $typ*, <$W x i1>)"
    push!(instrs, "call void @llvm.masked.compressstore.$(suffix(W,T))($vtyp %0, $typ* %ptr, <$W x i1> %mask)")
    push!(instrs, "ret void")
    quote
        $(Expr(:meta, :inline))
        Base.llvmcall(
            $((decl, join(instrs, "\n"))),
            Cvoid, Tuple{Vec{$W,$T}, Ptr{$T}, $U}, v, ptr, mask
        )
    end    
end
@inline compressstore!(ptr::Ptr{T}, v::Vec{W,T}, m::Mask{W}) where {W,T} = compressstore!(ptr, v, m.u)

@generated function expandload!(
    ::Type{Vec{W,T}}, ptr::Ptr{T}, mask::U
) where {W,T,U<:Unsigned}
    @assert 8sizeof(U) >= W
    ptyp = JuliaPointerType
    typ = llvmtype(T)
    vtyp = "<$W x $typ>"
    vptrtyp = "<$W x $typ*>"
    mtyp_input = llvmtype(U)
    mtyp_trunc = "i$W"
    instrs = String[]
    push!(instrs, "%ptr = inttoptr $ptyp %0 to $typ*")
    if mtyp_input == mtyp_trunc
        push!(instrs, "%mask = bitcast $mtyp_input %1 to <$W x i1>")
    else
        push!(instrs, "%masktrunc = trunc $mtyp_input %1 to $mtyp_trunc")
        push!(instrs, "%mask = bitcast $mtyp_trunc %masktrunc to <$W x i1>")
    end
    decl = "declare $vtyp @llvm.masked.expandload.$(suffix(W,T))($typ*, <$W x i1>, $vtyp)"
    push!(instrs, "%res = call $vtyp @llvm.masked.expandload.$(suffix(W,T))($typ* %ptr, <$W x i1> %mask, $vtyp zeroinitializer)")
    push!(instrs, "ret $vtyp %res")
    quote
        $(Expr(:meta, :inline))
        Base.llvmcall(
            $((decl, join(instrs, "\n"))),
            Vec{$W,$T}, Tuple{Ptr{$T}, $U}, ptr, mask
        )
    end
end
@inline expandload!(::Type{Vec{W,T}}, ptr::Ptr{T}, m::Mask{W}) where {W,T} = compressstore!(Vec{W,T}, ptr, m.u)


### LLVM syntax errors on invariants
# struct Invariant{L,T}
#     ivp::Ptr{Cvoid}
#     ptr::Ptr{T}
# end
# @generated function invariant_start!(ptr::Ptr{T}, ::Val{L}) where {L,T}
#     ptyp = JuliaPointerType
#     decls = "declare {}* @llvm.invariant.start.p0i8(i64, i8* nocapture)"
#     instrs = [
#         "%ptr = inttoptr $ptyp %0 to i8*",
#         "%ivt = call {}* @llvm.invariant.start.p0i8(i64 $(L*sizeof(T)), i8* %ptr)",
#         "%ivp = ptrtoint {}* %ivt to $ptyp",
#         "ret %ivp"
#     ]
#     quote
#         $(Expr(:meta,:inline))
#         ivp = Base.llvmcall(
#             $((decls, join(instrs, "\n"))),
#             Ptr{Cvoid}, Tuple{Ptr{$T}}, ptr
#         )
#         Invariant{$(L*sizeof(T)),T}(ivp, ptr)
#     end
# end
# @generated function invariant_end!(ivp::Invariant{L}) where {L}
#     ptyp = JuliaPointerType
#     decls = "declare void @llvm.invariant.end.p0i8({}*, i64, i8* nocapture)"
#     instrs = [
#         "%ivp = inttoptr $ptyp %0 to {}*",
#         "%ptr = inttoptr $ptyp %1 to i8*",
#         "call void @llvm.lifetime.end.p0i8({}* %ivp, i64 $(L), i8* %ptr)",
#         "ret void"
#     ]
#     quote
#         $(Expr(:meta,:inline))
#         Base.llvmcall(
#             $((decls, join(instrs, "\n"))),
#             Cvoid, Tuple{Ptr{$T}}, ivp.ivp, ivp.ptr
#         )
#     end
# end

for store ∈ [:vstore!, :vnoaliasstore!]
    @eval @inline $store(ptr::VectorizationBase.AbstractStridedPointer{T}, v::AbstractSIMDVector{W,T}, i, b::Bool) where {W,T} = (b && $store(ptr, v, i))

    @eval @inline $store(ptr::VectorizationBase.AbstractPointer{T1}, v::AbstractStructVec{W,T2}, i::Tuple) where {W,T1,T2} = $store(ptr, vconvert(Vec{W,T1}, v), i)
    @eval @inline $store(ptr::VectorizationBase.AbstractPointer{T1}, v::AbstractStructVec{W,T2}, i::Tuple, u::Unsigned) where {W,T1,T2} = $store(ptr, vconvert(Vec{W,T1}, v), i, u)
    @eval @inline $store(ptr::VectorizationBase.AbstractPointer{T1}, v::AbstractStructVec{W,T2}, i::Tuple, u::AbstractMask{W}) where {W,T1,T2} = $store(ptr, vconvert(Vec{W,T1}, v), i, tomask(u).u)

    @eval @inline $store(ptr::VectorizationBase.AbstractPointer{T}, m::AbstractMask{W}, i::Tuple) where {W,T} = $store(ptr, vifelse(tomask(m), vone(Vec{W,T}), vzero(Vec{W,T})), i)
    @eval @inline $store(ptr::VectorizationBase.AbstractPointer{T}, m::AbstractMask{W}, i::Tuple, mask::AbstractMask{W}) where {W,T} = $store(ptr, vifelse(tomask(m), vone(Vec{W,T}), vzero(Vec{W,T})), i, tomask(mask))
    @eval @inline function $store(ptr::VectorizationBase.AbstractPointer{Bool}, m::SVec{W,Bool}, i::Tuple) where {W}
        $store(ptr.ptr, extract_data(m), VectorizationBase.offset(ptr, VectorizationBase.staticm1(i)))
    end
    @eval @inline function $store(ptr::VectorizationBase.AbstractPointer{Bool}, m::SVec{W,Bool}, i::Tuple, mask::AbstractMask{W}) where {W}
        $store(ptr.ptr, extract_data(m), VectorizationBase.offset(ptr, VectorizationBase.staticm1(i)), tomask(mask))
    end
end
@inline vstore!(ptr::Ptr{T}, v::_MM{W}, i) where {W, T <: Integer} = vstore!(ptr, vrange(_MM{W}(v.i % T)), i)
@inline vstore!(ptr::Ptr{T}, v::_MM{W}, i, m::Mask) where {W, T <: Integer} = vstore!(ptr, vrange(_MM{W}(v.i % T)), i, m)


using VectorizationBase: AbstractColumnMajorStridedPointer, PackedStridedPointer, tdot
@inline VectorizationBase.gep(ptr::AbstractColumnMajorStridedPointer, i::NTuple{W,Core.VecElement{I}}) where {W,I<:Integer} = gep(ptr.ptr, i)

@inline vadd(s::I1, v::Vec{W,I2}) where {I1<:Integer,I2<:Integer,W} = vadd(vconvert(Vec{W,I2}, s), v)

# @inline VectorizationBase.gep(ptr::AbstractColumnMajorStridedPointer, i::Tuple) = @inbounds gep(ptr.ptr, vadd(i[1], tdot(ptr.strides, Base.tail(i))))
# @inline VectorizationBase.gep(ptr::AbstractColumnMajorStridedPointer{T,0}, i::Tuple) where {T} = @inbounds gep(ptr.ptr, i[1])
@inline VectorizationBase.tdot(a::Tuple{I,Any}, b::Tuple{Vec{W,I},Any}) where {W,I} = @inbounds vmul(first(a),first(b)) + tdot(Base.tail(a),Base.tail(b))
@inline VectorizationBase.tdot(a::Tuple{I,Any}, b::Tuple{SVec{W,I},Any}) where {W,I} = @inbounds vmul(first(a),extract_data(first(b))) + tdot(Base.tail(a),Base.tail(b))
@inline VectorizationBase.tdot(a::Tuple{I}, b::Tuple{Vec{W,I}}) where {W,I} = @inbounds vmul(first(a),first(b))
@inline VectorizationBase.tdot(a::Tuple{I}, b::Tuple{SVec{W,I}}) where {W,I} = @inbounds vmul(first(a),extract_data(first(b)))



@inline function VectorizationBase.tdot(a::Tuple{I1,Any}, b::Tuple{Vec{W,I2},Any}) where {W,I1,I2}
    @inbounds vmul(first(a) % I2,first(b)) + vconvert(Vec{W,I2},tdot(Base.tail(a),Base.tail(b)))
end
@inline function VectorizationBase.tdot(a::Tuple{I1,Any}, b::Tuple{SVec{W,I2},Any}) where {W,I1,I2}
    @inbounds vmul(first(a) % I2,extract_data(first(b))) + vconvert(Vec{W,I2},tdot(Base.tail(a),Base.tail(b)))
end
@inline function VectorizationBase.tdot(a::Tuple{I1}, b::Tuple{Vec{W,I2}}) where {W,I1,I2}
    @inbounds vmul(first(a) % I2,first(b))
end
@inline function VectorizationBase.tdot(a::Tuple{I1}, b::Tuple{SVec{W,I2}}) where {W,I1,I2}
    @inbounds vmul(first(a) % I2,extract_data(first(b)))
end

# _MM support
# zero initialized
# scalar if only and Int

for store ∈ [:vstore!, :vnoaliasstore!]
    @eval @inline $store(ptr::Ptr{T}, v::S, i::Union{SVec{W,<:Integer},_MM{W}}) where {W,T<:Number,S<:Number} = $store(ptr, vbroadcast(Vec{W,T}, v), i)
    @eval @inline $store(ptr::Ptr{T}, v::S, i::Union{SVec{W,<:Integer},_MM{W}}, u::Unsigned) where {W,T<:Number,S<:Number} = $store(ptr, vbroadcast(Vec{W,T}, v), i, u)
    @eval @inline $store(ptr::Ptr{T}, v::S, i::Union{SVec{W,<:Integer},_MM{W}}, u::Mask{W}) where {W,T<:Number,S<:Number} = $store(ptr, vbroadcast(Vec{W,T}, v), i, u.u)
end

vectypewidth(::Type{V}) where {W, V<:AbstractSIMDVector{W}} = W::Int
vectypewidth(::Type{_MM{W}}) where {W} = W::Int
vectypewidth(::Any) = 1

@inline vload(r::AbstractRange, i::Tuple{_MM{W,T}}) where {W,T} = SVec(vadd(vrangemul(Val{W}(), step(r), Val{0}()), @inbounds r[vadd(i[1].i, one(T))]))
@inline vload(r::UnitRange, i::Tuple{_MM{W,T}}) where {W,T} = @inbounds(_MM{W}(r[vadd(i[1].i, one(T))]))
# Ignore masks
@inline vload(r::AbstractRange, i::Tuple{_MM{W,T}}, ::Unsigned) where {W,T} = SVec(vadd(vrangemul(Val{W}(), step(r), Val{0}()), @inbounds r[vadd(i[1].i, one(T))]))
@inline vload(r::UnitRange, i::Tuple{_MM{W,T}}, ::Unsigned) where {W,T} = @inbounds(_MM{W}(r[vadd(i[1].i, one(T))]))
@inline vload(r::AbstractRange, i::Tuple{_MM{W,T}}, ::Mask) where {W,T} = SVec(vadd(vrangemul(Val{W}(), step(r), Val{0}()), @inbounds r[vadd(i[1].i, one(T))]))
@inline vload(r::UnitRange, i::Tuple{_MM{W,T}}, ::Mask) where {W,T} = @inbounds(_MM{W}(r[vadd(i[1].i, one(T))]))

function transposeshuffle0(split, W)
    tup = Expr(:tuple)
    w = 0
    S = 1 << split
    while w < W
        for s ∈ 0:S-1
            push!(tup.args, w + s)
        end
        for s ∈ 0:S-1
            # push!(tup.args, w + S + W + s)
            push!(tup.args, w + W + s)
        end
        w += 2S
    end
    Expr(:call, Expr(:curly, :Val, tup))
end
function transposeshuffle1(split, W)
    tup = Expr(:tuple)
    w = 0
    S = 1 << split
    while w < W
        for s ∈ 0:S-1
            push!(tup.args, w+S + s)
        end
        for s ∈ 0:S-1
            # push!(tup.args, w + W + s)
            push!(tup.args, w + S + W + s)
        end
        w += 2S
    end
    Expr(:call, Expr(:curly, :Val, tup))
end

@generated function vhaddstore!(ptr::AbstractPointer{T}, v::NTuple{N,V}, i::NTuple{D,I}) where {T,N,W,V<:AbstractSIMDVector{W,T},D,I<:Integer}
    q = Expr(:block, Expr(:meta, :inline), Expr(:(=), :bptr, Expr(:call, :gep, :ptr, :i)))
    if N > 1 && ispow2(N) && ispow2(W)
        extractblock = Expr(:block)
        vectors = [Symbol(:v_, n) for n ∈ 0:N-1]
        for n ∈ 1:N
            push!(extractblock.args, Expr(:(=), vectors[n], Expr(:ref, :v, n)))
        end
        push!(q.args, Expr(:macrocall, Symbol("@inbounds"), LineNumberNode(@__LINE__, Symbol(@__FILE__)), extractblock))
        ncomp = 0
        minWN = min(W,N)
        while ncomp < N
            Nt = minWN;
            Wt = W
            splits = 0
            while Nt > 1
                Nt >>>= 1
                shuffle0 = transposeshuffle0(splits, Wt)
                shuffle1 = transposeshuffle1(splits, Wt)
                splits += 1
                for nh ∈ 1:Nt
                    n1 = 2nh
                    n0 = n1 - 1
                    v0 = vectors[n0 + ncomp]; v1 = vectors[n1 + ncomp]; vh = vectors[nh + ncomp];
                    # combine n0 and n1
                    push!(q.args, Expr(
                        :(=), vh, Expr(
                            :call, :vadd,
                            Expr(:call, :shufflevector, v0, v1, shuffle0),
                            Expr(:call, :shufflevector, v0, v1, shuffle1))
                    )
                          )
                end
            end
            # v0 is now the only vector
            v0 = vectors[ncomp + 1]
            while Wt > minWN
                Wh = Wt >>> 1
                push!(q.args, Expr(
                    :(=), v0, Expr(
                        :call, :vadd,
                        Expr(:call, :shufflevector, v0, Expr(:call, Expr(:curly, :Val, Expr(:tuple, [w for w ∈ 0:Wh-1]...)))),
                        Expr(:call, :shufflevector, v0, Expr(:call, Expr(:curly, :Val, Expr(:tuple, [w for w ∈ Wh:Wt-1]...)))))
                )
                      )
                Wt = Wh
            end
            push!(q.args, Expr(:call, :vstore!, :bptr, v0))
            ncomp += minWN
        end
    else
        for n ∈ 1:N
            push!(
                q.args,
                Expr(
                    :call, :vstore!,
                    Expr(:call, :gep, :bptr, Expr(:tuple, n-1)),
                    Expr(
                        :call, :vsum,
                        Expr(:macrocall, Symbol("@inbounds"), LineNumberNode(@__LINE__, Symbol(@__FILE__)), Expr(:ref, :v, n))
                    )
                )
            )
        end
    end
    q
end

using VectorizationBase: MappedStridedPointer
@inline function vload(::Type{Vec{W,T}}, m::MappedStridedPointer{F,T}) where {W,F,T}
    extract_data(m.f(vload(SVec{W,T}, m.ptr)))
end
@inline function vload(::Type{Vec{W,T}}, m::MappedStridedPointer{F,T}, i) where {W,F,T}
    extract_data(m.f(vload(SVec{W,T}, m.ptr, i)))
end
@inline function vload(::Type{Vec{W,T}}, m::MappedStridedPointer{F,T}, mask::Union{Unsigned,Mask{W,<:Unsigned}}) where {W,F,T}
    extract_data(m.f(vload(SVec{W,T}, m.ptr, mask)))
end
@inline function vload(::Type{Vec{W,T}}, m::MappedStridedPointer{F,T}, i, mask::Union{Unsigned,Mask{W,<:Unsigned}}) where {W,F,T}
    extract_data(m.f(vload(SVec{W,T}, m.ptr, i, mask)))
end
@inline function vload(::Type{SVec{W,T}}, m::MappedStridedPointer{F,T}) where {W,F,T}
    m.f(vload(SVec{W,T}, m.ptr))
end
@inline function vload(::Type{SVec{W,T}}, m::MappedStridedPointer{F,T}, i) where {W,F,T}
    m.f(vload(SVec{W,T}, m.ptr, i))
end
@inline function vload(::Type{SVec{W,T}}, m::MappedStridedPointer{F,T}, mask::Union{Unsigned,Mask{W,<:Unsigned}}) where {W,F,T}
    m.f(vload(SVec{W,T}, m.ptr, mask))
end
@inline function vload(::Type{SVec{W,T}}, m::MappedStridedPointer{F,T}, i, mask::Union{Unsigned,Mask{W,<:Unsigned}}) where {W,F,T}
    m.f(vload(SVec{W,T}, m.ptr, i, mask))
end
@inline function vload(::Val{W}, m::MappedStridedPointer{F,T}) where {W,F,T}
    m.f(vload(SVec{W,T}, m.ptr))
end
@inline function vload(::Val{W}, m::MappedStridedPointer{F,T}, i) where {W,F,T}
    m.f(vload(SVec{W,T}, m.ptr, i))
end
@inline function vload(::Val{W}, m::MappedStridedPointer{F,T}, mask::Union{Unsigned,Mask{W,<:Unsigned}}) where {W,F,T}
    m.f(vload(SVec{W,T}, m.ptr, mask))
end
@inline function vload(::Val{W}, m::MappedStridedPointer{F,T}, i, mask::Union{Unsigned,Mask{W,<:Unsigned}}) where {W,F,T}
    m.f(vload(SVec{W,T}, m.ptr, i, mask))
end


@generated function vload(ptr::Ptr{Bool}, i::_MM{W,I}) where {W,I<:Integer}
    U = VectorizationBase.mask_type(W)
    utype = "i$(8sizeof(U))"
    itype = "i$(8sizeof(I))"
    ptyp = "i$(8sizeof(Int))"
    instrs = String[]
    push!(instrs, "%ptrbool = inttoptr $ptyp %0 to i8*")
    push!(instrs, "%ptroffset = getelementptr inbounds i8, i8* %ptrbool, $itype %1")
    push!(instrs, "%ptrvbool = bitcast i8* %ptroffset to <$W x i8>*")
    push!(instrs, "%vbool = load <$W x i8>, <$W x i8>* %ptrvbool")
    push!(instrs, "%vbooltrunc = trunc <$W x i8> %vbool to <$W x i1>")
    if 8sizeof(U) == W
        push!(instrs, "%mask = bitcast <$W x i1> %vbooltrunc to $utype")
    else
        push!(instrs, "%masktrunc = bitcast <$W x i1> %vbooltrunc to i$(W)")
        push!(instrs, "%mask = zext i$(W) %masktrunc to $utype")
    end
    push!(instrs, "ret $utype %mask")
    quote
        $(Expr(:meta,:inline))
        Mask{$W}(Base.llvmcall(
            $(join(instrs,"\n")),
            $U,
            Tuple{Ptr{Bool},$I},
            ptr, i.i
        ))
    end
end

@generated function vstore!(ptr::Ptr{Bool}, v::Mask{W,U}, i::I) where {W,U,I<:Integer}
    utype = "i$(8sizeof(U))"
    itype = "i$(8sizeof(I))"
    ptyp = "i$(8sizeof(Int))"
    instrs = String[]
    push!(instrs, "%ptrbool = inttoptr $ptyp %0 to i8*")
    push!(instrs, "%ptroffset = getelementptr inbounds i8, i8* %ptrbool, $itype %2")
    # push!(instrs, "%ptr = bitcast i8* %ptroffset to <$(8W) x i1>*")
    push!(instrs, "%ptr = bitcast i8* %ptroffset to <$(W) x i8>*")
    if W == 8sizeof(U)
        push!(instrs, "%mask = bitcast $utype %1 to <$W x i1>")
    else
        push!(instrs, "%masktrunc = trunc $utype %1 to i$(W)")
        push!(instrs, "%mask = bitcast i$(W) %masktrunc to <$W x i1>")
    end
    push!(instrs, "%maskzext = zext <$W x i1> %mask to <$W x i8>")
    # push!(instrs, "%vbits = bitcast <$W x i8> %maskzext to <$(8W) x i1>")
    # undefmask = ""
    # for w in 1:W
    #     for b in 1:7
    #         undefmask *= " i1 undef,"
    #     end
    #     undefmask *= w == W ? " i1 1" : " i1 1,"
    # end
    # push!(instrs, "%vundefbits = and <$(8W) x i1> %vbits, <$undefmask >")
    push!(instrs, "store <$(W) x i8> %maskzext, <$(W) x i8>* %ptr")
    push!(instrs, "ret void")
    quote
        $(Expr(:meta,:inline))
        Base.llvmcall(
            $(join(instrs,"\n")),
            Cvoid,
            Tuple{Ptr{Bool}, $U, $I},
            ptr, v.u, i
        )
    end
end




