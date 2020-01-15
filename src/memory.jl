
"""
I'm not sure on the details, but I think this function can only allocate up to

524288 doubles

that is it may not be able to allocate arrays with more than about half a million elements, or 4 megabytes.
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
@inline function alloca(N::Int, ::Type{T} = Float64, ::Val{Align} = Val{64}()) where {T, Align}
    alloca(Base.unsafe_trunc(Int32, N), T, Val{Align}())
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
    ::Union{Type{Vec{N,T}},Val{N}}, ptr::Ptr{T},
    ::Val{Aligned}, ::Val{Nontemporal}
    # ::Val{Aligned} = Val{false}(), ::Val{Nontemporal} = Val{false}()
) where {N,T,Aligned, Nontemporal}
    @assert isa(Aligned, Bool)
    ptyp = JuliaPointerType
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
    ::Union{Type{Vec{N,T}},Val{N}}, ptr::Ptr{T}, mask::Vec{N,Bool}, ::Val{Aligned}# = Val{false}()
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
    ::Union{Type{Vec{N,T}},Val{N}}, ptr::Ptr{T}, mask::U, ::Val{Aligned}# = Val{false}()
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

for v ∈ (:Vec, :SVec, :Val)
    vargs = Union{Symbol,Expr}[v === :Val ? :(::Val{W}) : :(::Type{$v{W,T}})]
    for ptr ∈ (:Ptr, :AbstractInitializedPointer)
        pargs = push!(copy(vargs), :(ptr::$ptr{T}))
        for index ∈ (true,false)
            if index
                icall = Union{Symbol,Expr}[:(Vec{W,T}), ptr == :Ptr ? :(ptr + i) : Expr(:call, :gep, :ptr, :i)]
                iargs = push!(copy(pargs), :i)
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
                    push!(fcall, Expr(:call, Expr(:curly, :Val, f !== :vload))) # aligned arg
                    if !mask
                        push!(fcall, Expr(:call, Expr(:curly, :Val, f === :vloadnt))) # nontemporal argf
                    end
                    body = Expr(:call, :vload, fcall...)
                    if v !== :Vec
                        body = Expr(:call, :SVec, body)
                    end
                    @eval @inline function $f($(margs...)) where {W,T}
                        $body
                    end
                end
            end
        end
    end
end

@inline vload(::Type{Vec{W,T}}, ::VectorizationBase.AbstractZeroInitializedPointer, args...) where {W,T} = vbroadcast(Vec{W,T}, zero(T))

@generated function vstore!(
    ptr::Ptr{T}, v::Vec{N,T},
    ::Val{Aligned}, ::Val{Nontemporal}
    # ::Val{Aligned} = Val{false}(), ::Val{Nontemporal} = Val{false}()
) where {N,T,Aligned, Nontemporal}
    @assert isa(Aligned, Bool)
    ptyp = JuliaPointerType
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
    ptyp = JuliaPointerType
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
@generated function vstore!(
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


for ptr ∈ (:Ptr, :AbstractPointer)#, :AbstractZeroInitializedPointer)
    pargs = Union{Symbol,Expr}[:(ptr::$ptr{T}), :(v::AbstractSIMDVector{W,T})]
    for index ∈ (true,false)
        if index
            icall = Union{Symbol,Expr}[ptr == :Ptr ? :(ptr + i) : Expr(:call, :gep, :ptr, :i), :(extract_data(v))]
            iargs = push!(copy(pargs), :i)
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
                push!(fcall, Expr(:call, Expr(:curly, :Val, f !== :vstore!))) # aligned arg
                if !mask
                    push!(fcall, Expr(:call, Expr(:curly, :Val, f === :vstorent!))) # nontemporal argf
                end
                body = Expr(:call, :vstore!, fcall...)
                @eval @inline function $f($(margs...)) where {W,T}
                    $body
                end
            end
        end
    end
end

# @generated function vloadscope(
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
#     ptyp = JuliaPointerType
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
    ptyp = JuliaPointerType
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
# @inline function vadd(ptr::Ptr{T}, inds::Vec{W,I}) where {W,T,I <: Intger}
#     gep(ptr, i)
# end
# @inline function vmuladd(s::I, inds::Vec{W,I}, ptr::Ptr{T}) where {W,T,I<:Union{UInt,Int}}
#     gep(
#     vreinterpret(Vec{W,Ptr{T}}, vmuladd(vbroadcast(Vec{W,I}, s), inds, vbroadcast(Vec{W,I}, reinterpret(I, ptr))))
# end

@inline scatter!(ptr::Ptr{T}, inds::AbstractSIMDVector{W,I}, v::AbstractSIMDVector{W,T}, mask::U) where {W,T,I<:IntegerTypes,U<:Unsigned} = scatter!(gep(ptr, extract_data(inds)), extract_data(v), mask)
@inline scatter!(ptr::Ptr{T}, inds::AbstractSIMDVector{W,I}, v::AbstractSIMDVector{W,T}) where {W,T,I<:IntegerTypes} = scatter!(gep(ptr, extract_data(inds)), extract_data(v))
@inline gather(ptr::Ptr{T}, inds::AbstractSIMDVector{W,I}, mask::U) where {W,T,I<:IntegerTypes,U<:Unsigned} = SVec(gather(gep(ptr,extract_data(inds)), mask))
@inline gather(ptr::Ptr{T}, inds::AbstractSIMDVector{W,I}) where {W,T,I<:IntegerTypes} = SVec(gather(gep(ptr,extract_data(inds))))
@inline gather(ptr::Ptr{T}, inds::Vec{W,I}, mask::U) where {W,T,I<:IntegerTypes,U<:Unsigned} = gather(gep(ptr,extract_data(inds)), mask)
@inline gather(ptr::Ptr{T}, inds::Vec{W,I}) where {W,T,I<:IntegerTypes} = gather(gep(ptr,extract_data(inds)))


@inline vload(::Type{Vec{W,T}}, A::AbstractArray, args...) where {W,T} = vload(Vec{W,T}, VectorizationBase.vectorizable(A), args...)
@inline gather(A::AbstractArray, args...) = gather(VectorizationBase.vectorizable(A), args...)
@inline vstore!(A::AbstractArray, args...) = vstore!(VectorizationBase.vectorizable(A), args...)
@inline scatter!(A::AbstractArray, args...) = scatter!(VectorizationBase.vectorizable(A), args...)

@inline function gather(
    ptr::VectorizationBase.AbstractInitializedStridedPointer{T},
    inds::AbstractSIMDVector{W,I}, mask::U
) where {T,W,I <: IntegerTypes,U}
    SVec(gather(gep(ptr, extract_data(inds)), mask))
end
@inline function gather(
    ptr::VectorizationBase.AbstractInitializedStridedPointer{T},
    inds::AbstractSIMDVector{W,I}
) where {T,W,I <: IntegerTypes}
    SVec(gather(gep(ptr, extract_data(inds))))
end
@inline function gather(
    ptr::VectorizationBase.AbstractInitializedStridedPointer{T},
    ci::Union{NTuple{N,Int},CartesianIndex{N}},
    inds::AbstractSIMDVector{W,I}, mask::U
) where {T,N,W,I <: IntegerTypes,U}
    SVec(gather(gep(gep(ptr,ci), extract_data(inds)), mask))
end
@inline function gather(
    ptr::VectorizationBase.AbstractInitializedStridedPointer{T},
    ci::Union{NTuple{N,Int},CartesianIndex{N}},
    inds::AbstractSIMDVector{W,I}
) where {T,N,W,I <: IntegerTypes}
    SVec(gather(gep(gep(ptr,ci),extract_data(inds))))
end
@inline function scatter!(
    ptr::VectorizationBase.AbstractStridedPointer{T},
    inds::AbstractSIMDVector{W,I}, v::Vec{W,T}, mask::U
) where {T,W,I <: IntegerTypes,U}
    scatter!(gep(ptr,extract_data(inds)), v, mask)
end
@inline function scatter!(
    ptr::VectorizationBase.AbstractStridedPointer{T},
    inds::AbstractSIMDVector{W,I}, v::Vec{W,T}
) where {T,W,I <: IntegerTypes}
    scatter!(gep(ptr,extract_data(inds)), v)
end
@inline function scatter!(
    ptr::VectorizationBase.AbstractStridedPointer{T},
    ci::Union{NTuple{N,Int},CartesianIndex{N}},
    inds::AbstractSIMDVector{W,I}, v::Vec{W,T}, mask::U
) where {T,N,W,I <: IntegerTypes,U}
    scatter!(gep(gep(ptr,ci),extract_data(inds)), v, mask)
end
@inline function scatter!(
    ptr::VectorizationBase.AbstractStridedPointer{T},
    ci::Union{NTuple{N,Int},CartesianIndex{N}},
    inds::AbstractSIMDVector{W,I}, v::Vec{W,T}
) where {T,N,W,I <: IntegerTypes}
    scatter!(gep(gep(ptr,ci),extract_data(inds)), v)
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

@inline function lifetime_start!(ptr::Pointer{T}) where {T}
    lifetime_start!(pointer(ptr), Val{1}())
end
@inline function lifetime_start!(ptr::Pointer{T}, ::Val{L}) where {T,L}
    lifetime_start!(pointer(ptr), Val{L}())
end
@inline function lifetime_start!(ptr::Ptr{T}) where {T}
    lifetime_start!(ptr, Val{1}())
end

@inline function lifetime_end!(ptr::Pointer{T}) where {T}
    lifetime_end!(pointer(ptr), Val{1}())
end
@inline function lifetime_end!(ptr::Pointer{T}, ::Val{L}) where {T,L}
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

function tuple_range_vector_expr(W)
    t = Expr(:tuple)
    for w ∈ zero(W):W-one(W)
        push!(t.args, Expr(:call, Expr(:(.), :Core, QuoteNode(:VecElement)), w))
    end
    t
end
function intrangetuple(W, ::Type{T}) where {T}
    if sizeof(T) == 8
        tuple_range_vector_expr(Int(W))
    elseif sizeof(T) == 4
        tuple_range_vector_expr(Int32(W))
    elseif sizeof(T) == 2
        tuple_range_vector_expr(Int16(W))
    elseif sizeof(T) == 1
        tuple_range_vector_expr(Int8(W))
    else
        tuple_range_vector_expr(Int(W))
    end
end
@generated function vrange(::Val{W}, ::Type{T}) where {W, T}
    Expr(:block, Expr(:meta,:inline), intrangetuple(W, T))
end
@generated function svrange(::Val{W}, ::Type{T}) where {W, T}
    Expr(:block, Expr(:meta,:inline), Expr(:call, :SVec, intrangetuple(W, T)))
end
@inline vrange(::Val{W}) where {W} = vrange(Val{W}(), Float64)
@inline svrange(::Val{W}) where {W} = svrange(Val{W}(), Float64)

@inline vload(::Type{Vec{W,T}}, ptr::VectorizationBase.AbstractInitializedStridedPointer, i) where {W,T} = vload(Vec{W,T}, gep(ptr, i))
@inline vload(::Type{Vec{W,T}}, ptr::VectorizationBase.AbstractInitializedStridedPointer, i, U::Unsigned) where {W,T} = vload(Vec{W,T}, gep(ptr, i), U)
@inline vstore!(ptr::VectorizationBase.AbstractStridedPointer{T}, v::T, i) where {T} = vstore!(gep(ptr, i), v)
@inline vstore!(ptr::VectorizationBase.AbstractStridedPointer{T}, v::T, i, U::Unsigned) where {T} = vstore!(gep(ptr, i), v, U)

@inline vstore!(ptr::VectorizationBase.AbstractPointer{T1}, v::AbstractSIMDVector{W,T2}, args...) where {W,T1,T2} = vstore!(ptr, vconvert(Vec{W,T1}, v), args...)
using VectorizationBase: AbstractPackedStridedPointer, PackedStridedPointer, tdot
@inline VectorizationBase.gep(ptr::AbstractPackedStridedPointer, i::NTuple{W,Core.VecElement{I}}) where {W,I<:Integer} = gep(ptr.ptr, i)

@inline vadd(s::I1, v::Vec{W,I2}) where {I1<:Integer,I2<:Integer,W} = vadd(vconvert(Vec{W,I2}, s), v)

@inline VectorizationBase.gep(ptr::AbstractPackedStridedPointer, i::Tuple) = @inbounds gep(ptr.ptr, vadd(first(i), tdot(ptr.strides, Base.tail(i))))
@inline VectorizationBase.tdot(a::Tuple{I,Any}, b::Tuple{Vec{W,I},Any}) where {W,I} = @inbounds vmul(first(a),first(b)) + tdot(Base.tail(a),Base.tail(b))
@inline VectorizationBase.tdot(a::Tuple{I,Any}, b::Tuple{SVec{W,I},Any}) where {W,I} = @inbounds vmul(first(a),extract_data(first(b))) + tdot(Base.tail(a),Base.tail(b))
@inline VectorizationBase.tdot(a::Tuple{I}, b::Tuple{Vec{W,I}}) where {W,I} = @inbounds vmul(first(a),first(b))
@inline VectorizationBase.tdot(a::Tuple{I}, b::Tuple{SVec{W,I}}) where {W,I} = @inbounds vmul(first(a),extract_data(first(b)))



@inline function VectorizationBase.tdot(a::Tuple{I1,Any}, b::Tuple{Vec{W,I2},Any}) where {W,I1,I2}
    @inbounds vmul(Base.unsafe_trunc(I2,first(a)),first(b)) + vconvert(Vec{W,I2},tdot(Base.tail(a),Base.tail(b)))
end
@inline function VectorizationBase.tdot(a::Tuple{I1,Any}, b::Tuple{SVec{W,I2},Any}) where {W,I1,I2}
    @inbounds vmul(Base.unsafe_trunc(I2,first(a)),extract_data(first(b))) + vconvert(Vec{W,I2},tdot(Base.tail(a),Base.tail(b)))
end
@inline function VectorizationBase.tdot(a::Tuple{I1}, b::Tuple{Vec{W,I2}}) where {W,I1,I2}
    @inbounds vmul(Base.unsafe_trunc(I2,first(a)),first(b))
end
@inline function VectorizationBase.tdot(a::Tuple{I1}, b::Tuple{SVec{W,I2}}) where {W,I1,I2}
    @inbounds vmul(Base.unsafe_trunc(I2,first(a)),extract_data(first(b)))
end


@inline function vload(::Val{W}, ptr::PackedStridedPointer{T}, i::Tuple{I1,Vec{W,I2}}) where {W,T,I1,I2}
    @inbounds SVec(gather(gep(ptr.ptr, vmuladd(Base.unsafe_trunc(I2,first(ptr.strides)), last(i), first(i))), Val{false}()))
end
@inline function vload(::Val{W}, ptr::PackedStridedPointer{T}, i::Tuple{I1,Vec{W,I2}}, mask::Unsigned) where {W,T,I1,I2}
    @inbounds SVec(gather(gep(ptr.ptr, vmuladd(Base.unsafe_trunc(I2,first(ptr.strides)), last(i), first(i))), mask, Val{false}()))
end
@inline function vstore!(ptr::AbstractPackedStridedPointer{T}, v::Vec{W,T}, i::Tuple{I1,Vec{W,I2}}) where {W,T,I1,I2}
    @inbounds scatter!(gep(ptr.ptr, vmuladd(Base.unsafe_trunc(I2,first(ptr.strides)), last(i), first(i))), v, Val{false}())
end
@inline function vstore!(ptr::AbstractPackedStridedPointer{T}, v::Vec{W,T}, i::Tuple{I1,Vec{W,I2}}, mask::Unsigned) where {W,T,I1,I2}
    @inbounds scatter!(gep(ptr.ptr, vmuladd(Base.unsafe_trunc(I2,first(ptr.strides)), last(i), first(i))), v, mask, Val{false}())
end

@inline vload(::Type{Vec{W,T}}, ptr::Vec{W,Ptr{T}}, ::Val{Aligned}, ::Val{Temporal}) where {W,T,Aligned,Temporal} = gather(ptr, Val{Aligned}())
@inline vload(::Type{Vec{W,T}}, ptr::Vec{W,Ptr{T}}, ::Val{Aligned}) where {W,T,Aligned} = gather(ptr, Val{Aligned}())
@inline vload(::Type{Vec{W,T}}, ptr::Vec{W,Ptr{T}}) where {W,T} = gather(ptr, Val{false}())
@inline vload(::Type{Vec{W,T}}, ptr::Vec{W,Ptr{T}}, mask::Unsigned, ::Val{Aligned}) where {W,T,Aligned} = gather(ptr, mask, Val{Aligned}())
@inline vload(::Type{Vec{W,T}}, ptr::Vec{W,Ptr{T}}, mask::Unsigned) where {W,T} = gather(ptr, mask, Val{false}())

@inline vload(::Val{W}, ptr::Vec{W,Ptr{T}}) where {W,T} = SVec(gather(ptr, Val{false}()))
@inline vload(::Val{W}, ptr::Vec{W,Ptr{T}}, mask::Unsigned) where {W,T} = SVec(gather(ptr, mask, Val{false}()))
# @inline vload(::Type{Vec{W,T}}, ptr::SVec{W,Ptr{T}}, ::Val{Aligned}, ::Val{Temporal}) where {W,T,Aligned,Temporal} = gather(ptr, Val{Aligned}())
# @inline vload(::Type{Vec{W,T}}, ptr::SVec{W,Ptr{T}}, ::Val{Aligned}) where {W,T,Aligned,Temporal} = gather(ptr, Val{Aligned}())
# @inline vload(::Type{Vec{W,T}}, ptr::SVec{W,Ptr{T}}) where {W,T,Aligned,Temporal} = gather(ptr, Val{false}())
# @inline vload(::Type{Vec{W,T}}, ptr::SVec{W,Ptr{T}}, mask::Unsigned, ::Val{Aligned}) where {W,T,Aligned} = gather(ptr, mask, Val{Aligned}())
# @inline vload(::Type{Vec{W,T}}, ptr::SVec{W,Ptr{T}}, mask::Unsigned) where {W,T} = gather(ptr, mask, Val{false}())

# @inline vload(::Type{Vec{W,T}}, ptr::AbstractSIMDVector{W,Ptr{T}}, ::Val{Aligned}, ::Val{Temporal}) where {W,T,Aligned,Temporal} = gather(ptr, Val{Aligned}())
# @inline vload(::Type{Vec{W,T}}, ptr::AbstractSIMDVector{W,Ptr{T}}, ::Val{Aligned}) where {W,T,Aligned,Temporal} = gather(ptr, Val{Aligned}())
# @inline vload(::Type{Vec{W,T}}, ptr::AbstractSIMDVector{W,Ptr{T}}) where {W,T,Aligned,Temporal} = gather(ptr, Val{false}())
# @inline vload(::Type{Vec{W,T}}, ptr::AbstractSIMDVector{W,Ptr{T}}, mask::Unsigned, ::Val{Aligned}) where {W,T,Aligned} = gather(ptr, mask, Val{Aligned}())
# @inline vload(::Type{Vec{W,T}}, ptr::AbstractSIMDVector{W,Ptr{T}}, mask::Unsigned) where {W,T} = gather(ptr, mask, Val{false}())

@inline vstore!(ptr::AbstractSIMDVector{W,Ptr{T}}, v::AbstractSIMDVector{W,T}, ::Val{Aligned}, ::Val{Temporal}) where {W,T,Aligned,Temporal} = scatter!(extract_data(ptr), extract_data(v), Val{Aligned}())
@inline vstore!(ptr::AbstractSIMDVector{W,Ptr{T}}, v::AbstractSIMDVector{W,T}, ::Val{Aligned}) where {W,T,Aligned} = scatter!(extract_data(ptr), extract_data(v), Val{Aligned}())
@inline vstore!(ptr::AbstractSIMDVector{W,Ptr{T}}, v::AbstractSIMDVector{W,T}) where {W,T} = scatter!(extract_data(ptr), extract_data(v), Val{false}())
@inline vstore!(ptr::AbstractSIMDVector{W,Ptr{T}}, v::AbstractSIMDVector{W,T}, mask::Unsigned, ::Val{Aligned}) where {W,T,Aligned} = scatter!(extract_data(ptr), extract_data(v), mask, Val{Aligned}())
@inline vstore!(ptr::AbstractSIMDVector{W,Ptr{T}}, v::AbstractSIMDVector{W,T}, mask::Unsigned) where {W,T} = scatter!(extract_data(ptr), extract_data(v), mask, Val{false}())


using VectorizationBase: stride1
for v ∈ (:Vec, :SVec, :Val)
    vargs = Union{Symbol,Expr}[v === :Val ? :(::Val{W}) : :(::Type{$v{W,T}}), :(ptr::VectorizationBase.SparseStridedPointer{T})]
    vcall = :(vmul(stride1(ptr), vrange(Val{W}())))
    for index ∈ (true,false)
        icall = copy(vcall)
        if index
            push!(icall.args, index ? Expr(:call, :gep, :ptr, :i) : :(ptr.ptr))
            iargs = index ? push!(copy(vargs), :i) : vargs
            for mask ∈ (true,false)
                mcall = Expr(:call, :gather, icall)
                margs = if mask
                    push!(mcall.args, :mask)
                    push!(copy(iargs), :(mask::Unsigned))
                else
                    iargs
                end
                push!(mcall.args, :(Val{false}()))
                body = v === :Vec ? mcall : Expr(:call, :SVec, mcall)
                @eval @inline function vload($(margs...)) where {W,T}
                    $body
                end
            end
        end
    end
end
@inline function vstore!(ptr::VectorizationBase.AbstractSparseStridedPointer{T}, v::AbstractSIMDVector{W,T}, i) where {W,T}
    scatter!(gep(gep(ptr, i), vmul(stride1(ptr), vrange(Val{W}()))), extract_data(v), Val{false}())
end
@inline function vstore!(ptr::VectorizationBase.AbstractSparseStridedPointer{T}, v::AbstractSIMDVector{W,T}, i, U::Unsigned) where {W,T}
    scatter!(gep(gep(ptr, i), vmul(stride1(ptr), vrange(Val{W}()))), extract_data(v), U, Val{false}())
end
@inline function vstore!(ptr::VectorizationBase.AbstractSparseStridedPointer{T}, v::AbstractSIMDVector{W,T}) where {W,T}
    scatter!(gep(ptr.ptr, vmul(stride1(ptr), vrange(Val{W}()))), extract_data(v), Val{false}())
end
@inline function vstore!(ptr::VectorizationBase.AbstractSparseStridedPointer{T}, v::AbstractSIMDVector{W,T}, U::Unsigned) where {W,T}
    scatter!(gep(ptr.ptr, vmul(stride1(ptr), vrange(Val{W}()))), extract_data(v), U, Val{false}())
end

using VectorizationBase: _MM, AbstractZeroInitializedPointer
@inline Base.:(+)(i::_MM{W}, j::AbstractSIMDVector{W}) where {W} = vadd(i.i, j)
@inline Base.:(+)(i::AbstractSIMDVector{W}, j::_MM{W}) where {W} = vadd(i, j.i)
@inline Base.:(*)(i::_MM{W}, j::AbstractSIMDVector{W}) where {W} = vmul(i.i, j)
@inline Base.:(*)(i::AbstractSIMDVector{W}, j::_MM{W}) where {W} = vmul(i, j.i)

# _MM support
# zero initialized
# scalar if only and Int
@inline vload(::AbstractZeroInitializedPointer{T}, ::Tuple{Int}) where {W,T} = zero(T)
# if a Vec of some kind, broadcast the zero
@inline vload(::AbstractZeroInitializedPointer{T}, ::Tuple{_MM{W},Vararg}) where {W,T} = vbroadcast(Val{W}(), zero(T))
@inline vload(::AbstractZeroInitializedPointer{T}, ::Tuple{V,Vararg}) where {W,T,I<:Integer,V<:AbstractSIMDVector{W,I}} = vbroadcast(Val{W}(), zero(T))
@inline vload(::AbstractZeroInitializedPointer{T}, ::Tuple{<:Integer,V,Vararg}) where {W,T,I<:Integer,V<:AbstractSIMDVector{W,I}} = vbroadcast(Val{W}, zero(T))
@inline vload(::AbstractZeroInitializedPointer{T}, ::Tuple{<:Integer,_MM{W},Vararg}) where {W,T} = vbroadcast(Val{W}, zero(T))
# peel off indices
@inline vload(ptr::AbstractZeroInitializedPointer{T}, i::Tuple{<:Integer,<:Integer,Vararg}) where {W,T} = vload(ptr, Base.tail(i))

vectypewidth(::Type{V}) where {W, V<:AbstractSIMDVector{W}} = W::Int
vectypewidth(::Type{_MM{W}}) where {W} = W::Int
vectypewidth(::Any) = 1
# returns expr for gep call, and bool-tuple (scalar,contiguous)
function packed_strided_ptr_index(Iparam, N, ::Type{T}) where {T}
    Ni = length(Iparam)
    @assert Ni == N+1
    Iₙ = Iparam[1]
    W::Int = vectypewidth(Iₙ)::Int
    indexpr = Expr(:ref, :i, 1)
    # check remaining indices.
    for n ∈ 2:Ni
        Iₙ = Iparam[n]
        Wₜ = vectypewidth(Iₙ)::Int
        W = W == 1 ? Wₜ : ((Wₜ == 1 || W == Wₜ) ? W : throw("$W ≠ $Wₜ but all vectors should be of the same width."))
        iexpr = Expr(:ref, :i, n)
        if Iₙ <: _MM
            iexpr = Expr(:call, :+, Expr(:call, :svrange, Expr(:call, Expr(:curly, :Val, W)), T), iexpr)
        end
        indexpr = Expr(:call, :muladd, iexpr, Expr(:ref, :s, n-1), indexpr)
    end
    indexpr = Expr(:macrocall, Symbol("@inbounds"), LineNumberNode(@__LINE__, @__FILE__), indexpr)
    W, Expr(:call, :gep, Expr(:(.), :ptr, QuoteNode(:ptr)), Expr(:call, :extract_data, indexpr))
end


@generated function vload(ptr::PackedStridedPointer{T,N}, i::I) where {T,N,I<:Tuple}
    W, gepcall = packed_strided_ptr_index(I.parameters, N, T)
    sexpr = Expr(:(=), :s, Expr(:(.), :ptr, QuoteNode(:strides)))
    if W == 1
        Expr(:block, Expr(:meta,:inline), sexpr, Expr(:call, :load, gepcall))
    else
        Expr(:block, Expr(:meta,:inline), sexpr, Expr(:call, :vload, Expr(:call, Expr(:curly, :Val, W)), gepcall))
    end
end
@generated function vload(ptr::PackedStridedPointer{T,N}, i::I, mask::Unsigned) where {T,N,I<:Tuple}
    W, gepcall = packed_strided_ptr_index(I.parameters, N, T)
    sexpr = Expr(:(=), :s, Expr(:(.), :ptr, QuoteNode(:strides)))
    if W == 1
        Expr(:block, Expr(:meta,:inline), sexpr, Expr(:call, :load, gepcall))
    else
        Expr(:block, Expr(:meta,:inline), sexpr, Expr(:call, :vload, Expr(:call, Expr(:curly, :Val, W)), gepcall, :mask))
    end
end
@generated function vstore!(ptr::AbstractPackedStridedPointer{T,N}, v::AbstractSIMDVector{W1,T}, i::I) where {W1,T,N,I<:Tuple}
    W, gepcall = packed_strided_ptr_index(I.parameters, N, T)
    sexpr = Expr(:(=), :s, Expr(:(.), :ptr, QuoteNode(:strides)))
    if W1 == 1
        Expr(:block, Expr(:meta,:inline), sexpr, Expr(:call, :store!, gepcall, Expr(:macrocall, Symbol("@inbounds"), LineNumberNode(@__LINE__, @__FILE__), Expr(:call, :firstval, :v))))
    else
        Expr(:block, Expr(:meta,:inline), sexpr, Expr(:call, :vstore!, gepcall, Expr(:call, :extract_data, :v)))
    end
end
@generated function vstore!(ptr::AbstractPackedStridedPointer{T,N}, v::AbstractSIMDVector{W1,T}, i::I, mask::Unsigned) where {W1,T,N,I<:Tuple}
    W, gepcall = packed_strided_ptr_index(I.parameters, N, T)
    sexpr = Expr(:(=), :s, Expr(:(.), :ptr, QuoteNode(:strides)))
    if W1 == 1
        Expr(:block, Expr(:meta,:inline), sexpr, Expr(:call, :store!, gepcall, Expr(:macrocall, Symbol("@inbounds"), LineNumberNode(@__LINE__, @__FILE__), Expr(:call, :firstval, :v))))
    else
        Expr(:block, Expr(:meta,:inline), sexpr, Expr(:call, :vstore!, gepcall, Expr(:call, :extract_data, :v), :mask))
    end
end

using VectorizationBase: SparseStridedPointer, AbstractSparseStridedPointer
@inline VectorizationBase.gep(ptr::AbstractSparseStridedPointer, i::NTuple{W,Core.VecElement{I}}) where {W,I<:Integer} = gep(ptr.ptr, vmul(i,vbroadcast(Val{W}(), first(ptr.strides))))
function sparse_strided_ptr_index(Iparam, N, ::Type{T}) where {T}
    Ni = length(Iparam)
    @assert Ni == N
    Iₙ = Iparam[1]
    W::Int = vectypewidth(Iₙ)::Int
    # indexpr = Expr(:ref, :i, 1)
    local indexpr::Expr
    # check remaining indices.
    for n ∈ 1:N
        Iₙ = Iparam[n]
        Wₜ = vectypewidth(Iₙ)::Int
        W = W == 1 ? Wₜ : ((Wₜ == 1 || W == Wₜ) ? W : throw("$W ≠ $Wₜ but all vectors should be of the same width."))
        iexpr = Expr(:ref, :i, n)
        if Iₙ <: _MM
            iexpr = Expr(:call, :+, Expr(:call, :svrange, Expr(:call, Expr(:curly, :Val, W)), T), iexpr)
        end
        indexpr = if n == 1
            Expr(:call, :*, iexpr, Expr(:ref, :s, n))
        else
            Expr(:call, :muladd, iexpr, Expr(:ref, :s, n), indexpr)
        end
    end
    indexpr = Expr(:macrocall, Symbol("@inbounds"), LineNumberNode(@__LINE__, @__FILE__), indexpr)
    W, Expr(:call, :gep, Expr(:(.), :ptr, QuoteNode(:ptr)), Expr(:call, :extract_data, indexpr))
end
@generated function vload(ptr::SparseStridedPointer{T,N}, i::I) where {T,N,I<:Tuple}
    W, gepcall = sparse_strided_ptr_index(I.parameters, N, T)
    sexpr = Expr(:(=), :s, Expr(:(.), :ptr, QuoteNode(:strides)))
    if W == 1
        Expr(:block, Expr(:meta,:inline), sexpr, Expr(:call, :load, gepcall))
    else
        Expr(:block, Expr(:meta,:inline), sexpr, Expr(:call, :vload, Expr(:call, Expr(:curly, :Val, W)), gepcall))
    end
end
@generated function vload(ptr::SparseStridedPointer{T,N}, i::I, mask::Unsigned) where {T,N,I<:Tuple}
    W, gepcall = sparse_strided_ptr_index(I.parameters, N, T)
    sexpr = Expr(:(=), :s, Expr(:(.), :ptr, QuoteNode(:strides)))
    if W == 1
        Expr(:block, Expr(:meta,:inline), sexpr, Expr(:call, :load, gepcall))
    else
        Expr(:block, Expr(:meta,:inline), sexpr, Expr(:call, :vload, Expr(:call, Expr(:curly, :Val, W)), gepcall, :mask))
    end
end
@generated function vstore!(ptr::AbstractSparseStridedPointer{T,N}, v::AbstractSIMDVector{W1,T}, i::I) where {W1,T,N,I<:Tuple}
    W, gepcall = sparse_strided_ptr_index(I.parameters, N, T)
    sexpr = Expr(:(=), :s, Expr(:(.), :ptr, QuoteNode(:strides)))
    if W1 == 1
        Expr(:block, Expr(:meta,:inline), sexpr, Expr(:call, :store!, gepcall, Expr(:macrocall, Symbol("@inbounds"), LineNumberNode(@__LINE__, @__FILE__), Expr(:call, :firstval, :v))))
    else
        Expr(:block, Expr(:meta,:inline), sexpr, Expr(:call, :vstore!, gepcall, Expr(:call, :extract_data, :v)))
    end
end
@generated function vstore!(ptr::AbstractSparseStridedPointer{T,N}, v::AbstractSIMDVector{W1,T}, i::I, mask::Unsigned) where {W1,T,N,I<:Tuple}
    W, gepcall = sparse_strided_ptr_index(I.parameters, N, T)
    sexpr = Expr(:(=), :s, Expr(:(.), :ptr, QuoteNode(:strides)))
    if W1 == 1
        Expr(:block, Expr(:meta,:inline), sexpr, Expr(:call, :store!, gepcall, Expr(:macrocall, Symbol("@inbounds"), LineNumberNode(@__LINE__, @__FILE__), Expr(:call, :firstval, :v))))
    else
        Expr(:block, Expr(:meta,:inline), sexpr, Expr(:call, :vstore!, gepcall, Expr(:call, :extract_data, :v), :mask))
    end
end

using VectorizationBase: StaticStridedPointer, AbstractStaticStridedPointer
@generated function VectorizationBase.gep(ptr::AbstractStaticStridedPointer{T,X}, i::NTuple{W,Core.VecElement{I}}) where {T,X,W,I<:Integer}
    s = first(X.parameters)::Int
    if s == 1
        Expr(:block, Expr(:meta, :inline), Expr(:call, :gep, Expr(:(.), :ptr, QuoteNode(:ptr)), :i))
    else
        Expr(:block, Expr(:meta, :inline), Expr(:call, :gep, Expr(:(.), :ptr, QuoteNode(:ptr)), Expr(:call, :vmul, Expr(:call,:vbroadcast, Expr(:call,Expr(:curly,:Val,W)), s), :i)))
    end
end
function static_strided_ptr_index(Iparam, Xparam, ::Type{T}) where {T}
    N = length(Iparam)
    @assert Ni == length(Xparam)
    Iₙ = Iparam[1]
    W::Int = vectypewidth(Iₙ)::Int
    # indexpr = Expr(:ref, :i, 1)
    local indexpr::Expr
    # check remaining indices.
    for n ∈ 1:N
        Iₙ = Iparam[n]
        Xₙ = (Xparam)::Int
        Wₜ = vectypewidth(Iₙ)::Int
        W = W == 1 ? Wₜ : ((Wₜ == 1 || W == Wₜ) ? W : throw("$W ≠ $Wₜ but all vectors should be of the same width."))
        iexpr = Expr(:ref, :i, n)
        if Iₙ <: _MM
            iexpr = Expr(:call, :+, Expr(:call, :svrange, Expr(:call, Expr(:curly, :Val, W)), T), iexpr)
        end
        indexpr = if n == 1
            Xₙ == 1 ? iexpr : Expr(:call, :*, iexpr, Xₙ)
        else
            if Xₙ == 1
                Expr(:call, :+, iexpr, indexpr)
            else
                Expr(:call, :muladd, iexpr, Xₙ, indexpr)
            end
        end
    end
    indexpr = Expr(:macrocall, Symbol("@inbounds"), LineNumberNode(@__LINE__, @__FILE__), indexpr)
    W, Expr(:call, :gep, Expr(:(.), :ptr, QuoteNode(:ptr)), Expr(:call, :extract_data, indexpr))
end
@generated function vload(ptr::StaticStridedPointer{T,X}, i::I) where {T,X,I<:Tuple}
    W, gepcall = static_strided_ptr_index(I.parameters, X.parameters, T)
    if W == 1
        Expr(:block, Expr(:meta,:inline), Expr(:call, :load, gepcall))
    else
        Expr(:block, Expr(:meta,:inline), Expr(:call, :vload, Expr(:call, Expr(:curly, :Val, W)), gepcall))
    end
end
@generated function vload(ptr::StaticStridedPointer{T,X}, i::I, mask::Unsigned) where {T,X,I<:Tuple}
    W, gepcall = static_strided_ptr_index(I.parameters, X.parameters, T)
    if W == 1
        Expr(:block, Expr(:meta,:inline), Expr(:call, :load, gepcall))
    else
        Expr(:block, Expr(:meta,:inline), Expr(:call, :vload, Expr(:call, Expr(:curly, :Val, W)), gepcall, :mask))
    end
end
@generated function vstore!(ptr::AbstractStaticStridedPointer{T,X}, v::AbstractSIMDVector{W1,T}, i::I) where {W1,T,X,I<:Tuple}
    W, gepcall = static_strided_ptr_index(I.parameters, X.parameters, T)
    if W1 == 1
        Expr(:block, Expr(:meta,:inline), Expr(:call, :store!, gepcall, Expr(:macrocall, Symbol("@inbounds"), LineNumberNode(@__LINE__, @__FILE__), Expr(:call, :firstval, :v))))
    else
        Expr(:block, Expr(:meta,:inline), Expr(:call, :vstore!, gepcall, Expr(:call, :extract_data, :v)))
    end
end
@generated function vstore!(ptr::AbstractStaticStridedPointer{T,X}, v::AbstractSIMDVector{W1,T}, i::I, mask::Unsigned) where {W1,T,X,I<:Tuple}
    W, gepcall = static_strided_ptr_index(I.parameters, X.parameters, T)
    if W1 == 1
        Expr(:block, Expr(:meta,:inline), Expr(:call, :store!, gepcall, Expr(:macrocall, Symbol("@inbounds"), LineNumberNode(@__LINE__, @__FILE__),  Expr(:call, :firstval, :v))))
    else
        Expr(:block, Expr(:meta,:inline), Expr(:call, :vstore!, gepcall, Expr(:call, :extract_data, :v), :mask))
    end
end

