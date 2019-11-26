
"""
I'm not sure on the details, but I think this function can only allocate up to

524288 doubles

that is it may not be able to allocate arrays with more than about half a million elements, or 4 megabytes.
"""
@generated function alloca(::Val{N}, ::Type{T} = Float64, ::Val{Align} = Val{64}()) where {N, T, Align}
    typ = llvmtype(T)
    ptyp = llvmtype(Int)
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

@generated function vload(
    ::Type{Vec{N,T}}, ptr::Ptr{T},
    ::Val{Aligned}, ::Val{Nontemporal}
    # ::Val{Aligned} = Val{false}(), ::Val{Nontemporal} = Val{false}()
) where {N,T,Aligned, Nontemporal}
    @assert isa(Aligned, Bool)
    ptyp = llvmtype(Int)
    typ = llvmtype(T)
    vtyp = "<$N x $typ>"
    decls = String[]
    instrs = String[]
    if Aligned
        align = N * sizeof(T)
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
@generated function vload(
    ::Type{Vec{N,T}}, ptr::Ptr{T}, mask::Vec{N,Bool}, ::Val{Aligned}# = Val{false}()
) where {N,T,Aligned}
    @assert isa(Aligned, Bool)
    ptyp = llvmtype(Int)
    typ = llvmtype(T)
    vtyp = "<$N x $typ>"
    btyp = llvmtype(Bool)
    vbtyp = "<$N x $btyp>"
    decls = String[]
    instrs = String[]
    if Aligned
        align = N * sizeof(T)
    else
        align = sizeof(T)   # This is overly optimistic
    end

    push!(instrs, "%ptr = inttoptr $ptyp %0 to $vtyp*")
    push!(instrs, "%mask = trunc $vbtyp %1 to <$N x i1>")
    push!(decls,
        "declare $vtyp @llvm.masked.load.$(suffix(N,T))($vtyp*, i32, <$N x i1>, $vtyp)"
    )
    push!(instrs,
        "%res = call $vtyp @llvm.masked.load.$(suffix(N,T))($vtyp* %ptr, i32 $align, <$N x i1> %mask, $vtyp $(llvmconst(N, T, 0)))"
    )
    push!(instrs, "ret $vtyp %res")
    quote
        $(Expr(:meta, :inline))
        Base.llvmcall($((join(decls, "\n"), join(instrs, "\n"))),
            Vec{$N,$T}, Tuple{Ptr{$T}, Vec{$N,Bool}}, ptr, mask)
    end
end
@generated function vload(
    ::Type{Vec{N,T}}, ptr::Ptr{T}, mask::U, ::Val{Aligned}# = Val{false}()
) where {N,T,Aligned,U<:Unsigned}
    @assert isa(Aligned, Bool)
    @assert 8sizeof(U) >= N
    ptyp = llvmtype(Int)
    typ = llvmtype(T)
    vtyp = "<$N x $typ>"
    mtyp_input = llvmtype(U)
    mtyp_trunc = "i$N"
    decls = String[]
    instrs = String[]
    if Aligned
        align = N * sizeof(T)
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
        "%res = call $vtyp @llvm.masked.load.$(suffix(N,T))($vtyp* %ptr, i32 $align, <$N x i1> %mask, $vtyp $(llvmconst(N, T, 0)))"
    )
    push!(instrs, "ret $vtyp %res")
    quote
        $(Expr(:meta, :inline))
        Base.llvmcall($((join(decls, "\n"), join(instrs, "\n"))),
            Vec{$N,$T}, Tuple{Ptr{$T}, $U}, ptr, mask)
    end
end

for v ∈ (:Vec, :SVec)
    vargs = [:(::Type{$v{W,T}})]
    for ptr ∈ (:Ptr, :Pointer, :ZeroInitializedPointer)
        pargs = push!(copy(vargs), :(ptr::$ptr{T}))
        for index ∈ (true,false)
            if index
                icall = Union{Symbol,Expr}[:(Vec{W,T}), ptr == :Ptr ? :(ptr + i) : :(ptr.ptr + sizeof(T)*i)]
                iargs = push!(copy(pargs), :(i::Int))
            else
                icall = Union{Symbol,Expr}[:(Vec{W,T}), ptr == :Ptr ? :ptr : :(ptr.ptr)]
                iargs = pargs
            end
            for mask ∈ (true,false)
                if mask
                    margs = push!(copy(iargs), :(mask::Union{Vec{W,Bool},<:Unsigned}))
                    mcall = push!(copy(icall), :mask)
                    fopts = (:vload,:vloada)
                else
                    margs = iargs
                    mcall = icall
                    fopts = (:vload,:vloada,:vloadnt)
                end
                for f ∈ fopts
                    fcall = copy(mcall)
                    if f == :vload # aligned arg
                        push!(fcall, :(Val{false}()))
                    else
                        push!(fcall, :(Val{true}()))
                    end
                    if !mask
                        if f == :vloadnt # nontemporal arg
                            push!(fcall, :(Val{true}()))
                        else
                            push!(fcall, :(Val{false}()))
                        end
                    end
                    if ptr === :ZeroInitializedPointer
                        body = Expr(:tuple, :vbroadcast, :(Vec{W,T}), :(zero(T)))
                    else
                        body = Expr(:call, :vload, fcall...)
                        if v == :SVec
                            body = :(SVec($body))
                        end
                    end
                    @eval @inline function $f($(margs...)) where {W,T}
                        $body
                    end
                end
            end
        end
    end
end

@generated function vstore!(
    ptr::Ptr{T}, v::Vec{N,T},
    ::Val{Aligned}, ::Val{Nontemporal}
    # ::Val{Aligned} = Val{false}(), ::Val{Nontemporal} = Val{false}()
) where {N,T,Aligned, Nontemporal}
    @assert isa(Aligned, Bool)
    ptyp = llvmtype(Int)
    typ = llvmtype(T)
    vtyp = "<$N x $typ>"
    decls = String[]
    instrs = String[]
    if Aligned# || Nontemporal
        align = N * sizeof(T)
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

@generated function vstore!(
    ptr::Ptr{T}, v::Vec{N,T}, mask::Vec{N,Bool}, ::Val{Aligned}# = Val{false}()
) where {N,T,Aligned}
    @assert isa(Aligned, Bool)
    ptyp = llvmtype(Int)
    typ = llvmtype(T)
    vtyp = "<$N x $typ>"
    btyp = llvmtype(Bool)
    vbtyp = "<$N x $btyp>"
    decls = String[]
    instrs = String[]
    if Aligned
        align = N * sizeof(T)
    else
        align = sizeof(T)   # This is overly optimistic
    end
    push!(instrs, "%ptr = inttoptr $ptyp %0 to $vtyp*")
    push!(instrs, "%mask = trunc $vbtyp %2 to <$N x i1>")
    push!(decls,
        "declare void @llvm.masked.store.$(suffix(N,T))($vtyp, $vtyp*, i32, <$N x i1>)")
    push!(instrs,
        "call void @llvm.masked.store.$(suffix(N,T))($vtyp %1, $vtyp* %ptr, i32 $align, <$N x i1> %mask)")
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
@generated function vstore!(
    ptr::Ptr{T}, v::Vec{N,T}, mask::U, ::Val{Aligned}# = Val{false}()
) where {N,T,Aligned,U<:Unsigned}
    @assert isa(Aligned, Bool)
    ptyp = llvmtype(Int)
    typ = llvmtype(T)
    vtyp = "<$N x $typ>"
    mtyp_input = llvmtype(U)
    mtyp_trunc = "i$N"
    decls = String[]
    instrs = String[]
    if Aligned
        align = N * sizeof(T)
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


for ptr ∈ (:Ptr, :Pointer, :ZeroInitializedPointer)
    pargs = [:(ptr::$ptr{T}), :(v::AbstractSIMDVector{W,T})]
    for index ∈ (true,false)
        if index
            icall = Union{Symbol,Expr}[ptr == :Ptr ? :(ptr + i) : :(ptr.ptr + sizeof(T)*i), :(extract_data(v))]
            iargs = push!(copy(pargs), :(i::Int))
        else
            icall = Union{Symbol,Expr}[ptr == :Ptr ? :ptr : :(ptr.ptr), :(extract_data(v))]
            iargs = pargs
        end
        for mask ∈ (true,false)
            if mask
                margs = push!(copy(iargs), :(mask::Union{Vec{W,Bool},<:Unsigned}))
                mcall = push!(copy(icall), :mask)
                fopts = (:vstore!,:vstorea!)
            else
                margs = iargs
                mcall = icall
                fopts = (:vstore!,:vstorea!,:vstorent!)
            end
            for f ∈ fopts
                fcall = copy(mcall)
                if f == :vstore! # aligned arg
                    push!(fcall, :(Val{false}()))
                else
                    push!(fcall, :(Val{true}()))
                end
                if !mask
                    if f == :vstorent! # nontemporal arg
                        push!(fcall, :(Val{true}()))
                    else
                        push!(fcall, :(Val{false}()))
                    end
                end
                body = Expr(:call, :vstore!, fcall...)
                @eval @inline function $f($(margs...)) where {W,T}
                    $body
                end
            end
        end
    end
end

@generated function vloadscope(
    ::Type{Vec{N,T}}, ptr::Ptr{T}, ::Val{Scope},
    ::Val{Aligned}, ::Val{Nontemporal}
) where {N, T, Scope, Aligned, Nontemporal}
    @assert isa(Aligned, Bool)
    ptyp = llvmtype(Int)
    typ = llvmtype(T)
    vtyp = "<$N x $typ>"
    domain,scope,list = Scope::NTuple{3,Int}
    # domain, scope, list = "\"domain\"", "\"scope\"", "\"list\""
    decls = String["!$domain = !{!$domain}","!$scope = !{!$scope, !$domain}", "!$list = !{!$scope}"]
    # decls = String[]
    instrs = String[]
    if Aligned
        align = N * sizeof(T)
    else
        align = sizeof(T)   # This is overly optimistic
    end
    flags = [""]
    align > 0 && push!(flags, "align $align")
    push!(flags, "!alias.scope !$list")
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
@generated function vstorescope!(
    ptr::Ptr{T}, v::Vec{N,T}, ::Val{Scope},
    ::Val{Aligned}, ::Val{Nontemporal}
    # ::Val{Aligned} = Val{false}(), ::Val{Nontemporal} = Val{false}()
) where {N,T,Scope,Aligned, Nontemporal}
    @assert isa(Aligned, Bool)
    ptyp = llvmtype(Int)
    typ = llvmtype(T)
    vtyp = "<$N x $typ>"
    domain,scope,list = Scope::NTuple{3,Int}
    decls = String["!$domain = !{!$domain}","!$scope = !{!$scope, !$domain}", "!$list = !{!$scope}"]
    # decls = String[]
    instrs = String[]
    if Aligned# || Nontemporal
        align = N * sizeof(T)
    else
        align = sizeof(T)   # This is overly optimistic
    end
    flags = [""]
    align > 0 && push!(flags, "align $align")
    push!(flags, "!noalias !$list")
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
function nonaliased_store_and_load!(storeptr::Ptr{T}, v::Vec{W,T}, loadptr::Ptr{T}) where {W,T}
    ptyp = llvmtype(Int)
    typ = llvmtype(T)
    vtyp = "<$W x $typ>"
    decls = """!0 = !{!0}
!1 = !{!1, !0}
!2 = !{!1}
"""
    instrs = """
%spt = inttoptr $ptyp %0 to $vtyp
%lpt = inttoptr $ptyp %2 to $vtyp
store $vtyp %1, $vtyp* %spt align $(sizeof(T)), !noalias !2
%res = load $vtyp, $vtyp* %lpt align $(sizeof(T)), !noalias !2
ret $vtyp %res
"""
    quote
        $(Expr(:meta, :inline))
        Base.llvmcall(
            ($decls, $instrs),
            Vec{$W,$T},
            Tuple{Ptr{$T},Vec{$W,$T},Ptr{$T}},
            storeptr, v, loadptr
        )
    end
end

@generated function vloadconstant(
    ::Type{Vec{N,T}}, ptr::Ptr{T},
    ::Val{Aligned}, ::Val{Nontemporal}
) where {N, T, Aligned, Nontemporal}
    @assert isa(Aligned, Bool)
    ptyp = llvmtype(Int)
    typ = llvmtype(T)
    vtyp = "<$N x $typ>"
    decls = String[]
    instrs = String[]
    if Aligned
        align = N * sizeof(T)
    else
        align = sizeof(T)   # This is overly optimistic
    end
    flags = [""]
    align > 0 && push!(flags, "align $align")
    push!(flags, "!tbaa !13")
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
# @generated function vloadscope(
#     ::Type{Vec{N,T}}, ptr::Ptr{T}, ::Val{Scope},
#     ::Val{Aligned}, ::Val{Nontemporal}
# ) where {N, T, Scope, Aligned, Nontemporal}
#     @assert isa(Aligned, Bool)
#     ptyp = llvmtype(Int)
#     typ = llvmtype(T)
#     vtyp = "<$N x $typ>"
#     domain,scope,list = Scope::NTuple{3,Int}
#     decls = String["!$domain = !{!$domain}","!$scope = !{!$scope, !$domain}", "!$list = !{!$scope}"]
#     # decls = String[]
#     instrs = String[]
#     if Aligned
#         align = N * sizeof(T)
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
# @generated function vstorescope!(
#     ptr::Ptr{T}, v::Vec{N,T}, ::Val{Scope},
#     ::Val{Aligned}, ::Val{Nontemporal}
#     # ::Val{Aligned} = Val{false}(), ::Val{Nontemporal} = Val{false}()
# ) where {N,T,Scope,Aligned, Nontemporal}
#     @assert isa(Aligned, Bool)
#     ptyp = llvmtype(Int)
#     typ = llvmtype(T)
#     vtyp = "<$N x $typ>"
#     domain,scope,list = Scope::NTuple{3,Int}
#     decls = String["!$domain = !{!$domain}","!$scope = !{!$scope, !$domain}", "!$list = !{!$scope}"]
#     # decls = String[]
#     instrs = String[]
#     if Aligned# || Nontemporal
#         align = N * sizeof(T)
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
    ptyp = llvmtype(Int)
    vptyp = "<$N x $ptyp>"
    typ = llvmtype(T)
    vtyp = "<$N x $typ>"
    vptrtyp = "<$N x $typ*>"
    decls = String[]
    instrs = String[]
    if Aligned
        align = N * sizeof(T)
    else
        align = sizeof(T)   # This is overly optimistic
    end
    mask = join((", i1 true" for i ∈ 2:N))
    push!(instrs, "%ptr = inttoptr $vptyp %0 to $vptrtyp")
    push!(decls,
        "declare $vtyp @llvm.masked.gather.$(suffix(N,T)).$(suffix(N,Ptr{T}))($vptrtyp, i32, <$N x i1>, $vtyp)")
    push!(instrs,
        "%res = call $vtyp @llvm.masked.gather.$(suffix(N,T)).$(suffix(N,Ptr{T}))($vptrtyp %ptr, i32 $align, <$N x i1> <i1 true$(mask)>, $vtyp undef)")#$(llvmconst(N, T, 0)))")
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
    ptyp = llvmtype(Int)
    vptyp = "<$N x $ptyp>"
    typ = llvmtype(T)
    vtyp = "<$N x $typ>"
    vptrtyp = "<$N x $typ*>"
    mtyp_input = llvmtype(U)
    mtyp_trunc = "i$N"
    decls = String[]
    instrs = String[]
    if Aligned
        align = N * sizeof(T)
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
    push!(decls,
        "declare $vtyp @llvm.masked.gather.$(suffix(N,T)).$(suffix(N,Ptr{T}))($vptrtyp, i32, <$N x i1>, $vtyp)")
    push!(instrs,
        "%res = call $vtyp @llvm.masked.gather.$(suffix(N,T)).$(suffix(N,Ptr{T}))($vptrtyp %ptr, i32 $align, <$N x i1> %mask, $vtyp $(llvmconst(N, T, 0)))")
    push!(instrs, "ret $vtyp %res")
    quote
        $(Expr(:meta, :inline))
        Base.llvmcall(
            $((join(decls, "\n"), join(instrs, "\n"))),
            Vec{$N,$T}, Tuple{Vec{$N,Ptr{$T}}, $U}, ptr, mask
        )
    end
end
@generated function scatter!(
    ptr::Vec{N,Ptr{T}}, v::Vec{N,T}, ::Val{Aligned} = Val{false}()
) where {N,T,Aligned}
    @assert isa(Aligned, Bool)
    ptyp = llvmtype(Int)
    vptyp = "<$N x $ptyp>"
    typ = llvmtype(T)
    vtyp = "<$N x $typ>"
    vptrtyp = "<$N x $typ*>"
    mtyp_input = llvmtype(U)
    mtyp_trunc = "i$N"
    decls = String[]
    instrs = String[]
    if Aligned
        align = N * sizeof(T)
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
            Cvoid, Tuple{Vec{$N,$T}, Vec{$N,Ptr{$T}}, $U}, v, ptr, mask)
    end
end
@generated function scatter!(
    ptr::Vec{N,Ptr{T}}, v::Vec{N,T}, mask::U, ::Val{Aligned} = Val{false}()
) where {N,T,Aligned,U<:Unsigned}
    @assert isa(Aligned, Bool)
    @assert 8sizeof(U) >= N
    ptyp = llvmtype(Int)
    vptyp = "<$N x $ptyp>"
    typ = llvmtype(T)
    vtyp = "<$N x $typ>"
    vptrtyp = "<$N x $typ*>"
    mtyp_input = llvmtype(U)
    mtyp_trunc = "i$N"
    decls = String[]
    instrs = String[]
    if Aligned
        align = N * sizeof(T)
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
@inline function vadd(ptr::Ptr{T}, inds::Vec{W,I}) where {W,T,I<:Union{UInt,Int}}
    pirate_reinterpret(Vec{W,Ptr{T}}, vadd(vbroadcast(Vec{W,I}, reinterpret(I, ptr)), inds))
end
@inline function vmuladd(s::I, inds::Vec{W,I}, ptr::Ptr{T}) where {W,T,I<:Union{UInt,Int}}
    pirate_reinterpret(Vec{W,Ptr{T}}, vmuladd(vbroadcast(Vec{W,I}, s), inds, vbroadcast(Vec{W,I}, reinterpret(I, ptr))))
end
@inline scatter!(ptr::Ptr{T}, inds::Vec{N,I}, v::Vec{N,T}, mask::U) where {N,T,I<:IntegerTypes,U<:Unsigned} = scatter!(vmuladd(sizeof(T), inds, ptr), v, mask)
@inline scatter!(ptr::Ptr{T}, inds::Vec{N,I}, v::Vec{N,T}) where {N,T,I<:IntegerTypes} = scatter!(vmuladd(sizeof(T), inds, ptr), v)
@inline gather(ptr::Ptr{T}, inds::Vec{N,I}, mask::U) where {N,T,I<:IntegerTypes,U<:Unsigned} = gather(vmuladd(sizeof(T),inds,ptr), mask)
@inline gather(ptr::Ptr{T}, inds::Vec{N,I}) where {N,T,I<:IntegerTypes} = gather(vmuladd(sizeof(T),inds,ptr))

@inline vload(::Type{Vec{W,T}}, A::AbstractArray, args...) where {W,T} = vload(Vec{W,T}, VectorizationBase.vectorizable(A), args...)
@inline gather(A::AbstractArray, args...) = gather(VectorizationBase.vectorizable(A), args...)
@inline vstore!(A::AbstractArray, args...) = vstore!(VectorizationBase.vectorizable(A), args...)
@inline scatter!(A::AbstractArray, args...) = scatter!(VectorizationBase.vectorizable(A), args...)

@generated function lifetime_start!(ptr::Ptr{T}, ::Val{L}) where {L,T}
    ptyp = llvmtype(Int)
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
    ptyp = llvmtype(Int)
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

@inline function lifetime_start!(ptr::Pointer{T}) where {T}
    SIMDPirates.lifetime_start!(pointer(ptr), Val{1}())
end
@inline function lifetime_start!(ptr::Ptr{T}) where {T}
    SIMDPirates.lifetime_start!(ptr, Val{1}())
end

@inline function lifetime_end!(ptr::Pointer{T}) where {T}
    SIMDPirates.lifetime_end!(pointer(ptr), Val{1}())
end
@inline function lifetime_end!(ptr::Ptr{T}) where {T}
    SIMDPirates.lifetime_end!(ptr, Val{1}())
end
# Fallback is to do nothing
@inline lifetime_start!(::Any) = nothing
@inline lifetime_end!(::Any) = nothing

@generated function noalias!(ptr::Ptr{T}) where {T}
    ptyp = llvmtype(Int)
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
    ptyp = llvmtype(Int)
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
    ptyp = llvmtype(Int)
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
    push!(instrs, "%res = call $vtyp @llvm.masked.expandload.$(suffix(N,T))($typ* %ptr, <$N x i1> %mask, $vtyp $(llvmconst(N, T, 0)))")
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
#     ptyp = llvmtype(Int)
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
#     ptyp = llvmtype(Int)
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




