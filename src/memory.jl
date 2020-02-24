


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
    decls = String[]
    instrs = String[]
    push!(instrs, "%ptr = alloca $typ, i32 $N, align $Align")
    push!(instrs, "%iptr = ptrtoint $typ* %ptr to $ptyp")
    # push!(instrs, "%ptr = alloca i8, i32 $(N*sizeof(T)), align $Align")
    # push!(instrs, "%iptr = ptrtoint i8* %ptr to $ptyp")
    push!(instrs, "ret $ptyp %iptr")
    quote
        $(Expr(:meta,:inline))
        Base.llvmcall(
            $((join(decls, "\n"), join(instrs, "\n"))),
            Ptr{$T}, Tuple{}
        )
    end
end
@generated function alloca(N::Int32, ::Type{T} = Float64, ::Val{Align} = Val{64}()) where {T, Align}
    typ = llvmtype(T)
    ptyp = JuliaPointerType
    decls = String[]
    instrs = String[]
    push!(instrs, "%ptr = alloca $typ, i32 %0, align $Align")
    push!(instrs, "%iptr = ptrtoint $typ* %ptr to $ptyp")
    # push!(instrs, "%ptr = alloca i8, i32 $(N*sizeof(T)), align $Align")
    # push!(instrs, "%iptr = ptrtoint i8* %ptr to $ptyp")
    push!(instrs, "ret $ptyp %iptr")
    quote
        $(Expr(:meta,:inline))
        Base.llvmcall(
            $((join(decls, "\n"), join(instrs, "\n"))),
            Ptr{$T}, Tuple{Int32}, N
        )
    end
end
@inline function alloca(N::Integer, ::Type{T} = Float64, ::Val{Align} = Val{64}()) where {T, Align}
    alloca(N % Int32, T, Val{Align}())
end


function valloc(::Type{T}, N::Int, sz::Int) where {T}
    @assert N > 0
    @assert sz >= 0
    # We use padding to align the address of the first element, and
    # also to ensure that we can access past the last element up to
    # the next full vector width
    padding = N-1 + mod(-sz, N)
    mem = Vector{T}(undef, sz + padding)
    addr = Int(pointer(mem))
    off = mod(-addr, N * sizeof(T))
    @assert mod(off, sizeof(T)) == 0
    off = fld(off, sizeof(T))
    @assert 0 <= off <= padding
    res = view(mem, off+1 : off+sz)
    addr2 = Int(pointer(res))
    @assert mod(addr2, N * sizeof(T)) == 0
    res
end
function valloc(f, ::Type{T}, N::Int, sz::Int) where {T}
    mem = valloc(T, N, sz)
    @inbounds for i in 1:sz
        mem[i] = f(i)
    end
    mem
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
@generated function load(
    ::Type{Vec{N,T}}, ptr::Ptr{T}, ::Val{Aligned}, ::Val{Nontemporal}
) where {N,T,Aligned, Nontemporal}
    @assert isa(Aligned, Bool)
    ptyp = JuliaPointerType
    typ = llvmtype(T)
    vtyp = "<$N x $typ>"
    decls = String[]
    instrs = String[]
    if Aligned
        align = Base.datatype_alignment(Vec{N,T})
    else
        align = sizeof(T)   # This is overly optimistic
    end
    flags = [""]
    align > 0 && push!(flags, "align $align")
    Nontemporal && push!(flags, "!nontemporal !{i32 1}")
    push!(instrs, "%ptr = inttoptr $ptyp %0 to $vtyp*")
    push!(instrs, "%res = load $vtyp, $vtyp* %ptr" * join(flags, ", "))
    push!(instrs, "ret $vtyp %res")
    quote
        $(Expr(:meta, :inline))
        Base.llvmcall($((join(decls, "\n"), join(instrs, "\n"))),
            Vec{$N,$T}, Tuple{Ptr{$T}}, ptr)
    end
end
@generated function load(
    ptr::Ptr{T}, i::_MM{W}, ::Val{Aligned}, ::Val{Nontemporal}
) where {W,T,Aligned, Nontemporal}
    @assert isa(Aligned, Bool)
    ptyp = JuliaPointerType
    typ = llvmtype(T)
    vtyp = "<$W x $typ>"
    decls = String[]
    instrs = String[]
    if Aligned
        align = Base.datatype_alignment(Vec{W,T})
    else
        align = sizeof(T)   # This is overly optimistic
    end
    flags = [""]
    align > 0 && push!(flags, "align $align")
    Nontemporal && push!(flags, "!nontemporal !{i32 1}")
    push!(instrs, "%typptr = inttoptr $ptyp %0 to $typ*")
    push!(instrs, "%offsetptr = getelementptr inbounds $typ, $typ* %typptr, $ptyp %1")
    push!(instrs, "%ptr = bitcast $typ* %offsetptr to $vtyp*")
    push!(instrs, "%res = load $vtyp, $vtyp* %ptr" * join(flags, ", "))
    push!(instrs, "ret $vtyp %res")
    quote
        $(Expr(:meta, :inline))
        SVec(Base.llvmcall($((join(decls, "\n"), join(instrs, "\n"))),
            Vec{$W,$T}, Tuple{Ptr{$T},Int}, ptr, i.i))
    end
end
# @inline function load(::Type{Vec{W,T}}, ptr::Ptr{T}, i::Integer, ::Val{Aligned}, ::Val{Nontemporal}) where {W,T,Aligned,Nontemporal}
    # extract_data(load(ptr, _MM{W}(i % Int), Val{Aligned}(), Val{Nontemporal}()))
# end
@generated function load(
    ::Type{Vec{N,T}}, ptr::Ptr{T}, mask::Vec{N,Bool}, ::Val{Aligned}
) where {N,T,Aligned}
    @assert isa(Aligned, Bool)
    ptyp = JuliaPointerType
    typ = llvmtype(T)
    vtyp = "<$N x $typ>"
    btyp = llvmtype(Bool)
    vbtyp = "<$N x $btyp>"
    decls = String[]
    instrs = String[]
    if Aligned
        align = Base.datatype_alignment(Vec{N,T})
    else
        align = sizeof(T)   # This is overly optimistic
    end

    push!(instrs, "%ptr = inttoptr $ptyp %0 to $vtyp*")
    push!(instrs, "%mask = trunc $vbtyp %1 to <$N x i1>")
    push!(decls,
        "declare $vtyp @llvm.masked.load.$(suffix(N,T))($vtyp*, i32, <$N x i1>, $vtyp)"
    )
    push!(instrs,
        "%res = call $vtyp @llvm.masked.load.$(suffix(N,T))($vtyp* %ptr, i32 $align, <$N x i1> %mask, $vtyp zeroinitializer)"#undef)"# 
    )
    push!(instrs, "ret $vtyp %res")
    quote
        $(Expr(:meta, :inline))
        Base.llvmcall($((join(decls, "\n"), join(instrs, "\n"))),
            Vec{$N,$T}, Tuple{Ptr{$T}, Vec{$N,Bool}}, ptr, mask)
    end
end
@generated function load(
    ::Type{Vec{N,T}}, ptr::Ptr{T}, mask::U, ::Val{Aligned}
) where {N,T,Aligned,U<:Unsigned}
    @assert isa(Aligned, Bool)
    @assert 8sizeof(U) >= N
    ptyp = JuliaPointerType
    typ = llvmtype(T)
    vtyp = "<$N x $typ>"
    mtyp_input = llvmtype(U)
    mtyp_trunc = "i$N"
    decls = String[]
    instrs = String[]
    if Aligned
        align = Base.datatype_alignment(Vec{N,T})
    else
        align = sizeof(T)   # This is overly optimistic
    end
    push!(instrs, "%ptr = inttoptr $ptyp %0 to $vtyp*")
    if mtyp_input == mtyp_trunc
        push!(instrs, "%mask = bitcast $mtyp_input %1 to <$N x i1>")
    else
        push!(instrs, "%masktrunc = trunc $mtyp_input %1 to $mtyp_trunc")
        push!(instrs, "%mask = bitcast $mtyp_trunc %masktrunc to <$N x i1>")
    end
    push!(decls,
        "declare $vtyp @llvm.masked.load.$(suffix(N,T))($vtyp*, i32, <$N x i1>, $vtyp)"
    )
    push!(instrs,
        "%res = call $vtyp @llvm.masked.load.$(suffix(N,T))($vtyp* %ptr, i32 $align, <$N x i1> %mask, $vtyp zeroinitializer)"#undef)"# 
    )
    push!(instrs, "ret $vtyp %res")
    quote
        $(Expr(:meta, :inline))
        Base.llvmcall($((join(decls, "\n"), join(instrs, "\n"))),
            Vec{$N,$T}, Tuple{Ptr{$T}, $U}, ptr, mask)
    end
end
@generated function load(
    ptr::Ptr{T}, i::_MM{W}, mask::U, ::Val{Aligned}
) where {W,T,U<:Unsigned,Aligned}
    ptyp = JuliaPointerType
    typ = llvmtype(T)
    vtyp = "<$W x $typ>"
    mtyp_input = llvmtype(U)
    mtyp_trunc = "i$W"
    decls = String[]
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
    push!(instrs,"%res = call $vtyp @llvm.masked.load.$(suffix(W,T))($vtyp* %ptr, i32 $align, <$W x i1> %mask, $vtyp zeroinitializer)")
    push!(instrs, "ret $vtyp %res")
    quote
        $(Expr(:meta, :inline))
        SVec(Base.llvmcall($((join(decls, "\n"), join(instrs, "\n"))),
            Vec{$W,$T}, Tuple{Ptr{$T}, Int, $U}, ptr, i.i, mask))
    end
end
# @inline function load(::Type{Vec{W,T}}, ptr::Ptr{T}, i::Integer, mask::Unsigned, ::Val{Aligned}) where {W,T,Aligned}
    # extract_data(load(ptr, _MM{W}(i % Int), mask, Val{Aligned}()))
# end

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
        for mask ∈ (true,false)
            if mask
                margs = push!(copy(iargs), :(mask::Unsigned))
                mcall = push!(copy(icall), :mask)
                fopts = [:load,:loada]
            else
                margs = iargs
                mcall = icall
                fopts = [:load,:loada,:loadnt]
            end
            for f ∈ fopts
                body = Expr(:call, :load, mcall...)
                push!(body.args, Expr(:call, Expr(:curly, :Val, f !== :load))) # aligned arg
                if !mask
                    push!(body.args, Expr(:call, Expr(:curly, :Val, f === :loadnt))) # nontemporal argf
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

@inline load(::Type{Vec{W,T}}, ::VectorizationBase.AbstractZeroInitializedPointer, args...) where {W,T} = vzero(Vec{W,T})

@generated function store!(
    ptr::Ptr{T}, v::Vec{N,T}, ::Val{Aligned}, ::Val{Nontemporal}
) where {N,T,Aligned, Nontemporal}
    @assert isa(Aligned, Bool)
    ptyp = JuliaPointerType
    typ = llvmtype(T)
    vtyp = "<$N x $typ>"
    decls = String[]
    instrs = String[]
    if Aligned# || Nontemporal
        align = Base.datatype_alignment(Vec{N,T})
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
            $((join(decls, "\n"), join(instrs, "\n"))),
            Cvoid, Tuple{Ptr{$T}, Vec{$N,$T}}, ptr, v
        )
    end
end
@generated function store!(
    ptr::Ptr{T}, v::Vec{N,T}, i::I, ::Val{Aligned}, ::Val{Nontemporal}
) where {N,T,I,Aligned, Nontemporal}
    @assert isa(Aligned, Bool)
    ityp = llvmtype(I)
    ptyp = JuliaPointerType
    typ = llvmtype(T)
    vtyp = "<$N x $typ>"
    decls = String[]
    instrs = String[]
    if Aligned# || Nontemporal
        align = Base.datatype_alignment(Vec{N,T})
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
            $((join(decls, "\n"), join(instrs, "\n"))),
            Cvoid, Tuple{Ptr{$T}, Vec{$N,$T}, $I}, ptr, v, i
        )
    end
end

@generated function store!(
    ptr::Ptr{T}, v::Vec{N,T}, mask::Vec{N,Bool}, ::Val{Aligned}# = Val{false}()
) where {N,T,Aligned}
    @assert isa(Aligned, Bool)
    ptyp = JuliaPointerType
    typ = llvmtype(T)
    vtyp = "<$N x $typ>"
    btyp = llvmtype(Bool)
    vbtyp = "<$N x $btyp>"
    decls = String[]
    instrs = String[]
    if Aligned
        align = Base.datatype_alignment(Vec{N,T})
    else
        align = sizeof(T)   # This is overly optimistic
    end
    push!(instrs, "%ptr = inttoptr $ptyp %0 to $vtyp*")
    push!(instrs, "%mask = trunc $vbtyp %2 to <$N x i1>")
    push!(decls, "declare void @llvm.masked.store.$(suffix(N,T))($vtyp, $vtyp*, i32, <$N x i1>)")
    push!(instrs, "call void @llvm.masked.store.$(suffix(N,T))($vtyp %1, $vtyp* %ptr, i32 $align, <$N x i1> %mask)")
    push!(instrs, "ret void")
    quote
        $(Expr(:meta, :inline))
        Base.llvmcall(
            $((join(decls, "\n"), join(instrs, "\n"))),
            Cvoid, Tuple{Ptr{$T}, Vec{$N,$T}, Vec{$N,Bool}},
            ptr, v, mask
        )
    end
end
@generated function store!(
    ptr::Ptr{T}, v::Vec{N,T}, mask::U, ::Val{Aligned}# = Val{false}()
) where {N,T,Aligned,U<:Unsigned}
    @assert isa(Aligned, Bool)
    ptyp = JuliaPointerType
    typ = llvmtype(T)
    vtyp = "<$N x $typ>"
    mtyp_input = llvmtype(U)
    mtyp_trunc = "i$N"
    decls = String[]
    instrs = String[]
    if Aligned
        align = Base.datatype_alignment(Vec{N,T})
    else
        align = sizeof(T)   # This is overly optimistic
    end
    push!(instrs, "%ptr = inttoptr $ptyp %0 to $vtyp*")
    if mtyp_input == mtyp_trunc
        push!(instrs, "%mask = bitcast $mtyp_input %2 to <$N x i1>")
    else
        push!(instrs, "%masktrunc = trunc $mtyp_input %2 to $mtyp_trunc")
        push!(instrs, "%mask = bitcast $mtyp_trunc %masktrunc to <$N x i1>")
    end
    push!(decls,
        "declare void @llvm.masked.store.$(suffix(N,T))($vtyp, $vtyp*, i32, <$N x i1>)"
    )
    push!(instrs,
        "call void @llvm.masked.store.$(suffix(N,T))($vtyp %1, $vtyp* %ptr, i32 $align, <$N x i1> %mask)"
    )
    push!(instrs, "ret void")
    quote
        $(Expr(:meta, :inline))
        Base.llvmcall($((join(decls, "\n"), join(instrs, "\n"))),
            Cvoid, Tuple{Ptr{$T}, Vec{$N,$T}, $U},
            ptr, v, mask)
    end
end
@generated function store!(
    ptr::Ptr{T}, v::Vec{N,T}, i::I, mask::U, ::Val{Aligned}# = Val{false}()
) where {N,T,Aligned,U<:Unsigned,I<:Integer}
    @assert isa(Aligned, Bool)
    ityp = llvmtype(I)
    ptyp = JuliaPointerType
    typ = llvmtype(T)
    vtyp = "<$N x $typ>"
    mtyp_input = llvmtype(U)
    mtyp_trunc = "i$N"
    decls = String[]
    instrs = String[]
    if Aligned
        align = Base.datatype_alignment(Vec{N,T})
    else
        align = sizeof(T)   # This is overly optimistic
    end
    push!(instrs, "%typptr = inttoptr $ptyp %0 to $typ*")
    push!(instrs, "%offsetptr = getelementptr inbounds $typ, $typ* %typptr, $ityp %2")
    push!(instrs, "%ptr = bitcast $typ* %offsetptr to $vtyp*")
    if mtyp_input == mtyp_trunc
        push!(instrs, "%mask = bitcast $mtyp_input %3 to <$N x i1>")
    else
        push!(instrs, "%masktrunc = trunc $mtyp_input %3 to $mtyp_trunc")
        push!(instrs, "%mask = bitcast $mtyp_trunc %masktrunc to <$N x i1>")
    end
    push!(decls,
        "declare void @llvm.masked.store.$(suffix(N,T))($vtyp, $vtyp*, i32, <$N x i1>)"
    )
    push!(instrs,
        "call void @llvm.masked.store.$(suffix(N,T))($vtyp %1, $vtyp* %ptr, i32 $align, <$N x i1> %mask)"
    )
    push!(instrs, "ret void")
    quote
        $(Expr(:meta, :inline))
        Base.llvmcall($((join(decls, "\n"), join(instrs, "\n"))),
            Cvoid, Tuple{Ptr{$T}, Vec{$N,$T}, $I, $U},
            ptr, v, i, mask)
    end
end


# for ptr ∈ (:Ptr, :AbstractPointer)#, :AbstractZeroInitializedPointer)
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
        for mask ∈ (true,false)
            if mask
                margs = push!(copy(iargs), :(mask::Union{Vec{W,Bool},<:Unsigned}))
                mcall = push!(copy(icall), :mask)
                fopts = (:store!,:storea!)
            else
                margs = iargs
                mcall = icall
                fopts = (:store!,:storea!,:storent!)
            end
            for f ∈ fopts
                body = Expr(:call, :store!)
                append!(body.args, mcall)
                push!(body.args, Expr(:call, Expr(:curly, :Val, f !== :store!))) # aligned arg
                if !mask
                    push!(body.args, Expr(:call, Expr(:curly, :Val, f === :storent!))) # nontemporal argf
                end
                @eval @inline function $f($(margs...)) where {W,T}
                    $body
                end
            end
        end
    end
end

# @generated function loadscope(
#     ::Type{Vec{N,T}}, ptr::Ptr{T}, ::Val{Scope},
#     ::Val{Aligned}, ::Val{Nontemporal}
# ) where {N, T, Scope, Aligned, Nontemporal}
#     @assert isa(Aligned, Bool)
#     ptyp = JuliaPointerType
#     typ = llvmtype(T)
#     vtyp = "<$N x $typ>"
#     domain,scope,list = Scope::NTuple{3,Int}
#     decls = String["!$domain = !{!$domain}","!$scope = !{!$scope, !$domain}", "!$list = !{!$scope}"]
#     # decls = String[]
#     instrs = String[]
#     if Aligned
#         align = Base.datatype_alignment(Vec{N,T})
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
#             Vec{$N,$T}, Tuple{Ptr{$T}}, ptr)
#     end
# end
# @generated function storescope!(
#     ptr::Ptr{T}, v::Vec{N,T}, ::Val{Scope},
#     ::Val{Aligned}, ::Val{Nontemporal}
#     # ::Val{Aligned} = Val{false}(), ::Val{Nontemporal} = Val{false}()
# ) where {N,T,Scope,Aligned, Nontemporal}
#     @assert isa(Aligned, Bool)
#     ptyp = JuliaPointerType
#     typ = llvmtype(T)
#     vtyp = "<$N x $typ>"
#     domain,scope,list = Scope::NTuple{3,Int}
#     decls = String["!$domain = !{!$domain}","!$scope = !{!$scope, !$domain}", "!$list = !{!$scope}"]
#     # decls = String[]
#     instrs = String[]
#     if Aligned# || Nontemporal
#         align = Base.datatype_alignment(Vec{N,T})
#     else
#         align = sizeof(T)   # This is overly optimistic
#     end
#     flags = [""]
#     align > 0 && push!(flags, "align $align")
#     push!(flags, "!noalias !$list")
#     Nontemporal && push!(flags, "!nontemporal !{i32 1}")
#     push!(instrs, "%ptr = inttoptr $ptyp %0 to $vtyp*")
#     push!(instrs, "store $vtyp %1, $vtyp* %ptr" * join(flags, ", "))
#     push!(instrs, "ret void")
#     quote
#         $(Expr(:meta, :inline))
#         Base.llvmcall(
#             $((join(decls, "\n"), join(instrs, "\n"))),
#             Cvoid, Tuple{Ptr{$T}, Vec{$N,$T}}, ptr, v
#         )
#     end
# end


@generated function gather(
   ptr::Vec{N,Ptr{T}}, ::Val{Aligned} = Val{false}()
) where {N,T,Aligned}
    @assert isa(Aligned, Bool)
    ptyp = JuliaPointerType
    vptyp = "<$N x $ptyp>"
    typ = llvmtype(T)
    vtyp = "<$N x $typ>"
    vptrtyp = "<$N x $typ*>"
    decls = String[]
    instrs = String[]
    if Aligned
        align = Base.datatype_alignment(Vec{N,T})
    else
        align = sizeof(T)   # This is overly optimistic
    end
    mask = join((", i1 true" for i ∈ 2:N))
    push!(instrs, "%ptr = inttoptr $vptyp %0 to $vptrtyp")
    push!(decls, "declare $vtyp @llvm.masked.gather.$(suffix(N,T)).$(suffix(N,Ptr{T}))($vptrtyp, i32, <$N x i1>, $vtyp)")
    push!(instrs, "%res = call $vtyp @llvm.masked.gather.$(suffix(N,T)).$(suffix(N,Ptr{T}))($vptrtyp %ptr, i32 $align, <$N x i1> <i1 true$(mask)>, $vtyp undef)")
    push!(instrs, "ret $vtyp %res")
    quote
        $(Expr(:meta, :inline))
        Base.llvmcall(
            $((join(decls, "\n"), join(instrs, "\n"))),
            Vec{$N,$T}, Tuple{Vec{$N,Ptr{$T}}}, ptr
        )
    end
end

@generated function gather(
   ptr::Vec{N,Ptr{T}}, mask::U, ::Val{Aligned} = Val{false}()
) where {N,T,Aligned,U<:Unsigned}
    @assert isa(Aligned, Bool)
    @assert 8sizeof(U) >= N
    ptyp = JuliaPointerType
    vptyp = "<$N x $ptyp>"
    typ = llvmtype(T)
    vtyp = "<$N x $typ>"
    vptrtyp = "<$N x $typ*>"
    mtyp_input = llvmtype(U)
    mtyp_trunc = "i$N"
    decls = String[]
    instrs = String[]
    if Aligned
        align = Base.datatype_alignment(Vec{N,T})
    else
        align = sizeof(T)   # This is overly optimistic
    end
    push!(instrs, "%ptr = inttoptr $vptyp %0 to $vptrtyp")
    if mtyp_input == mtyp_trunc
        push!(instrs, "%mask = bitcast $mtyp_input %1 to <$N x i1>")
    else
        push!(instrs, "%masktrunc = trunc $mtyp_input %1 to $mtyp_trunc")
        push!(instrs, "%mask = bitcast $mtyp_trunc %masktrunc to <$N x i1>")
    end
    push!(decls, "declare $vtyp @llvm.masked.gather.$(suffix(N,T)).$(suffix(N,Ptr{T}))($vptrtyp, i32, <$N x i1>, $vtyp)")
    push!(instrs, "%res = call $vtyp @llvm.masked.gather.$(suffix(N,T)).$(suffix(N,Ptr{T}))($vptrtyp %ptr, i32 $align, <$N x i1> %mask, $vtyp zeroinitializer)")#undef)")
    push!(instrs, "ret $vtyp %res")
    quote
        $(Expr(:meta, :inline))
        Base.llvmcall(
            $((join(decls, "\n"), join(instrs, "\n"))),
            Vec{$N,$T}, Tuple{Vec{$N,Ptr{$T}}, $U}, ptr, mask
        )
    end
end
@generated function load(
   ptr::Ptr{T}, i::Vec{N,I}, ::Val{Aligned} = Val{false}()
) where {N,T,I<:Integer,Aligned}
    @assert isa(Aligned, Bool)
    ptyp = JuliaPointerType
    vptyp = "<$N x $ptyp>"
    typ = llvmtype(T)
    vtyp = "<$N x $typ>"
    vptrtyp = "<$N x $typ*>"
    decls = String[]
    instrs = String[]
    if Aligned
        align = Base.datatype_alignment(Vec{N,T})
    else
        align = sizeof(T)   # This is overly optimistic
    end
    ityp = llvmtype(I)
    vityp = "<$N x $ityp>"
    push!(instrs, "%sptr = inttoptr $ptyp %0 to $typ*")
    push!(instrs, "%ptr = getelementptr inbounds $typ, $typ* %sptr, $vityp %1")
    mask = join((", i1 true" for i ∈ 2:N))
    push!(decls, "declare $vtyp @llvm.masked.gather.$(suffix(N,T)).$(suffix(N,Ptr{T}))($vptrtyp, i32, <$N x i1>, $vtyp)")
    push!(instrs, "%res = call $vtyp @llvm.masked.gather.$(suffix(N,T)).$(suffix(N,Ptr{T}))($vptrtyp %ptr, i32 $align, <$N x i1> <i1 true$(mask)>, $vtyp undef)")
    push!(instrs, "ret $vtyp %res")
    quote
        $(Expr(:meta, :inline))
        Base.llvmcall(
            $((join(decls, "\n"), join(instrs, "\n"))),
            Vec{$N,$T}, Tuple{Ptr{$T},Vec{$N,$I}}, ptr, i
        )
    end
end
@inline load(ptr::Ptr, i::SVec{W,I}, ::Val{Aligned}) where {W,I<:Integer,Aligned} = SVec(load(ptr, extract_data(i), Val{Aligned}()))
@inline load(ptr::Ptr, i::SVec{W,I}, mask::Unsigned, ::Val{Aligned}) where {W,I<:Integer,Aligned} = SVec(load(ptr, extract_data(i), mask, Val{Aligned}()))
@inline load(ptr::Ptr, i::SVec{W,I}, ::Val{Aligned}, ::Val{false}) where {W,I<:Integer,Aligned} = SVec(load(ptr, extract_data(i), Val{Aligned}()))
@inline load(ptr::Ptr, i::SVec{W,I}, mask::Unsigned, ::Val{Aligned}, ::Val{false}) where {W,I<:Integer,Aligned} = SVec(load(ptr, extract_data(i), mask, Val{Aligned}()))
@generated function load(
   ptr::Ptr{T}, i::Vec{N,I}, mask::U, ::Val{Aligned} = Val{false}()
) where {N,T,I<:Integer,Aligned,U<:Unsigned}
    @assert isa(Aligned, Bool)
    @assert 8sizeof(U) >= N
    ptyp = JuliaPointerType
    vptyp = "<$N x $ptyp>"
    typ = llvmtype(T)
    vtyp = "<$N x $typ>"
    vptrtyp = "<$N x $typ*>"
    mtyp_input = llvmtype(U)
    mtyp_trunc = "i$N"
    decls = String[]
    instrs = String[]
    if Aligned
        align = Base.datatype_alignment(Vec{N,T})
    else
        align = sizeof(T)   # This is overly optimistic
    end
    ityp = llvmtype(I)
    vityp = "<$N x $ityp>"
    push!(instrs, "%sptr = inttoptr $ptyp %0 to $typ*")
    push!(instrs, "%ptr = getelementptr inbounds $typ, $typ* %sptr, $vityp %1")
    if mtyp_input == mtyp_trunc
        push!(instrs, "%mask = bitcast $mtyp_input %2 to <$N x i1>")
    else
        push!(instrs, "%masktrunc = trunc $mtyp_input %2 to $mtyp_trunc")
        push!(instrs, "%mask = bitcast $mtyp_trunc %masktrunc to <$N x i1>")
    end
    push!(decls, "declare $vtyp @llvm.masked.gather.$(suffix(N,T)).$(suffix(N,Ptr{T}))($vptrtyp, i32, <$N x i1>, $vtyp)")
    push!(instrs, "%res = call $vtyp @llvm.masked.gather.$(suffix(N,T)).$(suffix(N,Ptr{T}))($vptrtyp %ptr, i32 $align, <$N x i1> %mask, $vtyp zeroinitializer)")#undef)")
    push!(instrs, "ret $vtyp %res")
    quote
        $(Expr(:meta, :inline))
        Base.llvmcall(
            $((join(decls, "\n"), join(instrs, "\n"))),
            Vec{$N,$T}, Tuple{Ptr{$T},Vec{$N,$I},$U}, ptr, i, mask
        )
    end
end
@generated function scatter!(
    ptr::Vec{N,Ptr{T}}, v::Vec{N,T}, ::Val{Aligned} = Val{false}()
) where {N,T,Aligned}
    @assert isa(Aligned, Bool)
    ptyp = JuliaPointerType
    vptyp = "<$N x $ptyp>"
    typ = llvmtype(T)
    vtyp = "<$N x $typ>"
    vptrtyp = "<$N x $typ*>"
    decls = String[]
    instrs = String[]
    if Aligned
        align = Base.datatype_alignment(Vec{N,T})
    else
        align = sizeof(T)   # This is overly optimistic
    end
    push!(instrs, "%ptr = inttoptr $vptyp %1 to $vptrtyp")
    mask = join((", i1 true" for i ∈ 2:N))
    # push!(decls, "declare void @llvm.masked.scatter.$(suffix(N,T)).$(suffix(N,Ptr{T}))($vtyp, $vptrtyp, i32, <$N x i1>)")
    # push!(instrs, "call void @llvm.masked.scatter.$(suffix(N,T)).$(suffix(N,Ptr{T}))($vtyp %0, $vptrtyp %ptr, i32 $align, <$N x i1> <i1 true$(mask)>)")
    push!(decls, "declare void @llvm.masked.scatter.$(suffix(N,T))($vtyp, $vptrtyp, i32, <$N x i1>)")
    push!(instrs, "call void @llvm.masked.scatter.$(suffix(N,T))($vtyp %0, $vptrtyp %ptr, i32 $align, <$N x i1> <i1 true$(mask)>)")
    push!(instrs, "ret void")
    quote
        $(Expr(:meta, :inline))
        Base.llvmcall($((join(decls, "\n"), join(instrs, "\n"))),
            Cvoid, Tuple{Vec{$N,$T}, Vec{$N,Ptr{$T}}}, v, ptr)
    end
end
@generated function scatter!(
    ptr::Vec{N,Ptr{T}}, v::Vec{N,T}, mask::U, ::Val{Aligned} = Val{false}()
) where {N,T,Aligned,U<:Unsigned}
    @assert isa(Aligned, Bool)
    @assert 8sizeof(U) >= N
    ptyp = JuliaPointerType
    vptyp = "<$N x $ptyp>"
    typ = llvmtype(T)
    vtyp = "<$N x $typ>"
    vptrtyp = "<$N x $typ*>"
    mtyp_input = llvmtype(U)
    mtyp_trunc = "i$N"
    decls = String[]
    instrs = String[]
    if Aligned
        align = Base.datatype_alignment(Vec{N,T})
    else
        align = sizeof(T)   # This is overly optimistic
    end
    push!(instrs, "%ptr = inttoptr $vptyp %1 to $vptrtyp")
    if mtyp_input == mtyp_trunc
        push!(instrs, "%mask = bitcast $mtyp_input %2 to <$N x i1>")
    else
        push!(instrs, "%masktrunc = trunc $mtyp_input %2 to $mtyp_trunc")
        push!(instrs, "%mask = bitcast $mtyp_trunc %masktrunc to <$N x i1>")
    end
    # push!(decls,
        # "declare void @llvm.masked.scatter.$(suffix(N,T)).$(suffix(N,Ptr{T}))($vtyp, $vptrtyp, i32, <$N x i1>)")
    # push!(instrs,
        # "call void @llvm.masked.scatter.$(suffix(N,T)).$(suffix(N,Ptr{T}))($vtyp %0, $vptrtyp %ptr, i32 $align, <$N x i1> %mask)")
    push!(decls, "declare void @llvm.masked.scatter.$(suffix(N,T))($vtyp, $vptrtyp, i32, <$N x i1>)")
    push!(instrs, "call void @llvm.masked.scatter.$(suffix(N,T))($vtyp %0, $vptrtyp %ptr, i32 $align, <$N x i1> %mask)")
    push!(instrs, "ret void")
    quote
        $(Expr(:meta, :inline))
        Base.llvmcall(
            $((join(decls, "\n"), join(instrs, "\n"))),
            Cvoid, Tuple{Vec{$N,$T}, Vec{$N,Ptr{$T}}, $U}, v, ptr, mask
        )
    end
end
@generated function store!(
    ptr::Ptr{T}, v::Vec{N,T}, i::Vec{N,I}, ::Val{Aligned} = Val{false}()
) where {N,T,Aligned, I<:Integer}
    @assert isa(Aligned, Bool)
    ptyp = JuliaPointerType
    vptyp = "<$N x $ptyp>"
    typ = llvmtype(T)
    vtyp = "<$N x $typ>"
    vptrtyp = "<$N x $typ*>"
    decls = String[]
    instrs = String[]
    if Aligned
        align = Base.datatype_alignment(Vec{N,T})
    else
        align = sizeof(T)   # This is overly optimistic
    end
    ityp = llvmtype(I)
    vityp = "<$N x $ityp>"
    push!(instrs, "%sptr = inttoptr $ptyp %0 to $typ*")
    push!(instrs, "%ptr = getelementptr inbounds $typ, $typ* %sptr, $vityp %2")
    mask = join((", i1 true" for i ∈ 2:N))
    # push!(decls, "declare void @llvm.masked.scatter.$(suffix(N,T)).$(suffix(N,Ptr{T}))($vtyp, $vptrtyp, i32, <$N x i1>)")
    # push!(instrs, "call void @llvm.masked.scatter.$(suffix(N,T)).$(suffix(N,Ptr{T}))($vtyp %0, $vptrtyp %ptr, i32 $align, <$N x i1> <i1 true$(mask)>)")
    push!(decls, "declare void @llvm.masked.scatter.$(suffix(N,T))($vtyp, $vptrtyp, i32, <$N x i1>)")
    push!(instrs, "call void @llvm.masked.scatter.$(suffix(N,T))($vtyp %1, $vptrtyp %ptr, i32 $align, <$N x i1> <i1 true$(mask)>)")
    push!(instrs, "ret void")
    quote
        $(Expr(:meta, :inline))
        Base.llvmcall($((join(decls, "\n"), join(instrs, "\n"))),
            Cvoid, Tuple{Ptr{$T}, Vec{$N,$T}, Vec{$N,$I}}, ptr, v, i)
    end
end
@generated function store!(
    ptr::Ptr{T}, v::Vec{N,T}, i::Vec{N,I}, mask::U, ::Val{Aligned} = Val{false}()
) where {N,T,Aligned,U<:Unsigned,I<:Integer}
    @assert isa(Aligned, Bool)
    @assert 8sizeof(U) >= N
    ptyp = JuliaPointerType
    vptyp = "<$N x $ptyp>"
    typ = llvmtype(T)
    vtyp = "<$N x $typ>"
    vptrtyp = "<$N x $typ*>"
    mtyp_input = llvmtype(U)
    mtyp_trunc = "i$N"
    decls = String[]
    instrs = String[]
    if Aligned
        align = Base.datatype_alignment(Vec{N,T})
    else
        align = sizeof(T)   # This is overly optimistic
    end
    ityp = llvmtype(I)
    vityp = "<$N x $ityp>"
    push!(instrs, "%sptr = inttoptr $ptyp %0 to $typ*")
    push!(instrs, "%ptr = getelementptr inbounds $typ, $typ* %sptr, $vityp %2")
    if mtyp_input == mtyp_trunc
        push!(instrs, "%mask = bitcast $mtyp_input %3 to <$N x i1>")
    else
        push!(instrs, "%masktrunc = trunc $mtyp_input %3 to $mtyp_trunc")
        push!(instrs, "%mask = bitcast $mtyp_trunc %masktrunc to <$N x i1>")
    end
    # push!(decls,
        # "declare void @llvm.masked.scatter.$(suffix(N,T)).$(suffix(N,Ptr{T}))($vtyp, $vptrtyp, i32, <$N x i1>)")
    # push!(instrs,
        # "call void @llvm.masked.scatter.$(suffix(N,T)).$(suffix(N,Ptr{T}))($vtyp %0, $vptrtyp %ptr, i32 $align, <$N x i1> %mask)")
    push!(decls, "declare void @llvm.masked.scatter.$(suffix(N,T))($vtyp, $vptrtyp, i32, <$N x i1>)")
    push!(instrs, "call void @llvm.masked.scatter.$(suffix(N,T))($vtyp %1, $vptrtyp %ptr, i32 $align, <$N x i1> %mask)")
    push!(instrs, "ret void")
    quote
        $(Expr(:meta, :inline))
        Base.llvmcall(
            $((join(decls, "\n"), join(instrs, "\n"))),
            Cvoid, Tuple{Ptr{$T}, Vec{$N,$T}, Vec{$N,$I}, $U}, ptr, v, i, mask
        )
    end
end

@generated function lifetime_start!(ptr::Ptr{T}, ::Val{L}) where {L,T}
    ptyp = JuliaPointerType
    decls = "declare void @llvm.lifetime.start(i64, i8* nocapture)"
    instrs = [
        "%ptr = inttoptr $ptyp %0 to i8*",
        "call void @llvm.lifetime.start(i64 $(L*sizeof(T)), i8* %ptr)",
        "ret void"
    ]
    quote
        $(Expr(:meta,:inline))
        Base.llvmcall(
            $((decls, join(instrs, "\n"))),
            Cvoid, Tuple{Ptr{$T}}, ptr
        )
    end
end
@generated function lifetime_end!(ptr::Ptr{T}, ::Val{L}) where {L,T}
    ptyp = JuliaPointerType
    decls = "declare void @llvm.lifetime.end(i64, i8* nocapture)"
    instrs = [
        "%ptr = inttoptr $ptyp %0 to i8*",
        "call void @llvm.lifetime.end(i64 $(L*sizeof(T)), i8* %ptr)",
        "ret void"
    ]
    quote
        $(Expr(:meta,:inline))
        Base.llvmcall(
            $((decls, join(instrs, "\n"))),
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

@generated function noalias!(ptr::Ptr{T}) where {T}
    ptyp = JuliaPointerType
    typ = llvmtype(T)
    decls = "define noalias $typ* @noalias($typ *%a) noinline { ret $typ* %a }"
    instrs = [
        "%ptr = inttoptr $ptyp %0 to $typ*",
        "%naptr = call $typ* @noalias($typ* %ptr)",
        "%jptr = ptrtoint $typ* %naptr to $ptyp",
        "ret $ptyp %jptr"
    ]
    quote
        $(Expr(:meta,:inline))
        Base.llvmcall(
            $((decls, join(instrs, "\n"))),
            Ptr{$T}, Tuple{Ptr{$T}}, ptr
        )
    end    
end

@generated function compressstore!(
    ptr::Ptr{T}, v::Vec{N,T}, mask::U
) where {N,T,U<:Unsigned}
    @assert 8sizeof(U) >= N
    ptyp = JuliaPointerType
    typ = llvmtype(T)
    vtyp = "<$N x $typ>"
    mtyp_input = llvmtype(U)
    mtyp_trunc = "i$N"
    decls = String[]
    instrs = String[]
    push!(instrs, "%ptr = inttoptr $ptyp %1 to $typ*")
    if mtyp_input == mtyp_trunc
        push!(instrs, "%mask = bitcast $mtyp_input %2 to <$N x i1>")
    else
        push!(instrs, "%masktrunc = trunc $mtyp_input %2 to $mtyp_trunc")
        push!(instrs, "%mask = bitcast $mtyp_trunc %masktrunc to <$N x i1>")
    end
    push!(decls, "declare void @llvm.masked.compressstore.$(suffix(N,T))($vtyp, $typ*, <$N x i1>)")
    push!(instrs, "call void @llvm.masked.compressstore.$(suffix(N,T))($vtyp %0, $typ* %ptr, <$N x i1> %mask)")
    push!(instrs, "ret void")
    quote
        $(Expr(:meta, :inline))
        Base.llvmcall(
            $((join(decls, "\n"), join(instrs, "\n"))),
            Cvoid, Tuple{Vec{$N,$T}, Ptr{$T}, $U}, v, ptr, mask
        )
    end    
end

@generated function expandload!(
    ::Type{Vec{N,T}}, ptr::Ptr{T}, mask::U
) where {N,T,U<:Unsigned}
    @assert 8sizeof(U) >= N
    ptyp = JuliaPointerType
    typ = llvmtype(T)
    vtyp = "<$N x $typ>"
    vptrtyp = "<$N x $typ*>"
    mtyp_input = llvmtype(U)
    mtyp_trunc = "i$N"
    decls = String[]
    instrs = String[]
    push!(instrs, "%ptr = inttoptr $ptyp %0 to $typ*")
    if mtyp_input == mtyp_trunc
        push!(instrs, "%mask = bitcast $mtyp_input %1 to <$N x i1>")
    else
        push!(instrs, "%masktrunc = trunc $mtyp_input %1 to $mtyp_trunc")
        push!(instrs, "%mask = bitcast $mtyp_trunc %masktrunc to <$N x i1>")
    end
    push!(decls, "declare $vtyp @llvm.masked.expandload.$(suffix(N,T))($typ*, <$N x i1>, $vtyp)")
    push!(instrs, "%res = call $vtyp @llvm.masked.expandload.$(suffix(N,T))($typ* %ptr, <$N x i1> %mask, $vtyp zeroinitializer)")
    push!(instrs, "ret $vtyp %res")
    quote
        $(Expr(:meta, :inline))
        Base.llvmcall(
            $((join(decls, "\n"), join(instrs, "\n"))),
            Vec{$N,$T}, Tuple{Ptr{$T}, $U}, ptr, mask
        )
    end
end


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


@inline store!(ptr::VectorizationBase.AbstractStridedPointer{T}, v::AbstractSIMDVector{W,T}, i, b::Bool) where {W,T} = (b && store!(ptr, v, i))

@inline store!(ptr::VectorizationBase.AbstractPointer{T1}, v::AbstractStructVec{W,T2}, i::Tuple) where {W,T1,T2} = store!(ptr, vconvert(Vec{W,T1}, v), i)
@inline store!(ptr::VectorizationBase.AbstractPointer{T1}, v::AbstractStructVec{W,T2}, i::Tuple, u::Unsigned) where {W,T1,T2} = store!(ptr, vconvert(Vec{W,T1}, v), i, u)
using VectorizationBase: AbstractPackedStridedPointer, PackedStridedPointer, tdot
@inline VectorizationBase.gep(ptr::AbstractPackedStridedPointer, i::NTuple{W,Core.VecElement{I}}) where {W,I<:Integer} = gep(ptr.ptr, i)

@inline vadd(s::I1, v::Vec{W,I2}) where {I1<:Integer,I2<:Integer,W} = vadd(vconvert(Vec{W,I2}, s), v)

@inline VectorizationBase.gep(ptr::AbstractPackedStridedPointer, i::Tuple) = @inbounds gep(ptr.ptr, vadd(first(i), tdot(ptr.strides, Base.tail(i))))
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
@inline load(::AbstractZeroInitializedPointer{T}, ::Tuple{Int}) where {W,T<:Number} = zero(T)
# if a Vec of some kind, broadcast the zero
@inline load(::AbstractZeroInitializedPointer{T}, ::Tuple{_MM{W},Vararg}) where {W,T<:Number} = vzero(SVec{W,T})
@inline load(::AbstractZeroInitializedPointer{T}, ::Tuple{V,Vararg}) where {W,T<:Number,I<:Integer,V<:AbstractSIMDVector{W,I}} = vzero(SVec{W,T})
@inline load(::AbstractZeroInitializedPointer{T}, ::Tuple{<:Integer,V,Vararg}) where {W,T<:Number,I<:Integer,V<:AbstractSIMDVector{W,I}} = vzero(SVec{W,T})
@inline load(::AbstractZeroInitializedPointer{T}, ::Tuple{<:Integer,_MM{W},Vararg}) where {W,T<:Number} = vzero(SVec{W,T})
# peel off indices
@inline load(ptr::AbstractZeroInitializedPointer{T}, i::Tuple{<:Integer,<:Integer,Vararg}) where {W,T<:Number} = load(ptr, Base.tail(i))

@inline store!(ptr::Ptr{T}, v::T, i::Union{SVec{W,<:Integer},_MM{W}}) where {W,T<:Number,I} = store!(ptr, vbroadcast(Vec{W,T}, v), i)
@inline store!(ptr::Ptr{T}, v::T, i::Union{SVec{W,<:Integer},_MM{W}}, u::Unsigned) where {W,T<:Number,I} = store!(ptr, vbroadcast(Vec{W,T}, v), i, u)

vectypewidth(::Type{V}) where {W, V<:AbstractSIMDVector{W}} = W::Int
vectypewidth(::Type{_MM{W}}) where {W} = W::Int
vectypewidth(::Any) = 1

@inline load(r::AbstractRange{T}, i::Tuple{_MM{W}}) where {W,T} = SVec(vadd(vrangemul(Val{W}(), step(r), Val{0}()), @inbounds r[i[1].i + 1]))
@inline load(r::UnitRange{T}, i::Tuple{_MM{W}}) where {W,T} = @inbounds(_MM{W}(r[i[1].i + 1]))
# Ignore masks
@inline load(r::AbstractRange{T}, i::Tuple{_MM{W}}, ::Unsigned) where {W,T} = SVec(vadd(vrangemul(Val{W}(), step(r), Val{0}()), @inbounds r[i[1].i + 1]))
@inline load(r::UnitRange{T}, i::Tuple{_MM{W}}, ::Unsigned) where {W,T} = @inbounds(_MM{W}(r[i[1].i + 1]))

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
            push!(q.args, Expr(:call, :store!, :bptr, v0))
            ncomp += minWN
        end
    else
        for n ∈ 1:N
            push!(
                q.args,
                Expr(
                    :call, :store!,
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
@inline function load(::Type{Vec{W,T}}, m::MappedStridedPointer{F,T}) where {W,F,T}
    extract_data(m.f(load(SVec{W,T}, m.ptr)))
end
@inline function load(::Type{Vec{W,T}}, m::MappedStridedPointer{F,T}, i) where {W,F,T}
    extract_data(m.f(load(SVec{W,T}, m.ptr, i)))
end
@inline function load(::Type{Vec{W,T}}, m::MappedStridedPointer{F,T}, mask::Union{<:Unsigned,Vec{W,Bool}}) where {W,F,T}
    extract_data(m.f(load(SVec{W,T}, m.ptr, mask)))
end
@inline function load(::Type{Vec{W,T}}, m::MappedStridedPointer{F,T}, i, mask::Union{<:Unsigned,Vec{W,Bool}}) where {W,F,T}
    extract_data(m.f(load(SVec{W,T}, m.ptr, i, mask)))
end
@inline function load(::Type{SVec{W,T}}, m::MappedStridedPointer{F,T}) where {W,F,T}
    m.f(load(SVec{W,T}, m.ptr))
end
@inline function load(::Type{SVec{W,T}}, m::MappedStridedPointer{F,T}, i) where {W,F,T}
    m.f(load(SVec{W,T}, m.ptr, i))
end
@inline function load(::Type{SVec{W,T}}, m::MappedStridedPointer{F,T}, mask::Union{<:Unsigned,Vec{W,Bool}}) where {W,F,T}
    m.f(load(SVec{W,T}, m.ptr, mask))
end
@inline function load(::Type{SVec{W,T}}, m::MappedStridedPointer{F,T}, i, mask::Union{<:Unsigned,Vec{W,Bool}}) where {W,F,T}
    m.f(load(SVec{W,T}, m.ptr, i, mask))
end
@inline function load(::Val{W}, m::MappedStridedPointer{F,T}) where {W,F,T}
    m.f(load(SVec{W,T}, m.ptr))
end
@inline function load(::Val{W}, m::MappedStridedPointer{F,T}, i) where {W,F,T}
    m.f(load(SVec{W,T}, m.ptr, i))
end
@inline function load(::Val{W}, m::MappedStridedPointer{F,T}, mask::Union{<:Unsigned,Vec{W,Bool}}) where {W,F,T}
    m.f(load(SVec{W,T}, m.ptr, mask))
end
@inline function load(::Val{W}, m::MappedStridedPointer{F,T}, i, mask::Union{<:Unsigned,Vec{W,Bool}}) where {W,F,T}
    m.f(load(SVec{W,T}, m.ptr, i, mask))
end



