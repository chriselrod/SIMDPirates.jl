


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
        "%res = call $vtyp @llvm.masked.load.$(suffix(N,T))($vtyp* %ptr, i32 $align, <$N x i1> %mask, $vtyp zeroinitializer)"#undef)"# 
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
        "%res = call $vtyp @llvm.masked.load.$(suffix(N,T))($vtyp* %ptr, i32 $align, <$N x i1> %mask, $vtyp zeroinitializer)"#undef)"# 
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

@inline vload(::Type{Vec{W,T}}, ::VectorizationBase.AbstractZeroInitializedPointer, args...) where {W,T} = vzero(Vec{W,T})

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


@inline vload(::Type{Vec{W,T}}, ptr::VectorizationBase.AbstractInitializedStridedPointer, i) where {W,T} = vload(Vec{W,T}, gep(ptr, i))
@inline vload(::Type{Vec{W,T}}, ptr::VectorizationBase.AbstractInitializedStridedPointer, i, U::Unsigned) where {W,T} = vload(Vec{W,T}, gep(ptr, i), U)
@inline vstore!(ptr::VectorizationBase.AbstractStridedPointer{T}, v::Vec{W,T}, i) where {W,T} = vstore!(gep(ptr, i), v)
@inline vstore!(ptr::VectorizationBase.AbstractStridedPointer{T}, v::Vec{W,T}, i, U::Unsigned) where {W,T} = vstore!(gep(ptr, i), v, U)
@inline vstore!(ptr::VectorizationBase.AbstractStridedPointer{T}, v::AbstractSIMDVector{W,T}, i, b::Bool) where {W,T} = (b && vstore!(ptr, v, i))

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


@inline function vload(::Val{W}, ptr::PackedStridedPointer{T}, i::Tuple{I1,Vec{W,I2}}) where {W,T,I1,I2}
    @inbounds SVec(gather(gep(ptr.ptr, vmuladd(first(ptr.strides) % I2, last(i), first(i))), Val{false}()))
end
@inline function vload(::Val{W}, ptr::PackedStridedPointer{T}, i::Tuple{I1,Vec{W,I2}}, mask::Unsigned) where {W,T,I1,I2}
    @inbounds SVec(gather(gep(ptr.ptr, vmuladd(first(ptr.strides) % I2, last(i), first(i))), mask, Val{false}()))
end
@inline function vstore!(ptr::AbstractPackedStridedPointer{T}, v::Vec{W,T}, i::Tuple{I1,Vec{W,I2}}) where {W,T,I1,I2}
    @inbounds scatter!(gep(ptr.ptr, vmuladd(first(ptr.strides) % I2, last(i), first(i))), v, Val{false}())
end
@inline function vstore!(ptr::AbstractPackedStridedPointer{T}, v::Vec{W,T}, i::Tuple{I1,Vec{W,I2}}, mask::Unsigned) where {W,T,I1,I2}
    @inbounds scatter!(gep(ptr.ptr, vmuladd(first(ptr.strides) % I2, last(i), first(i))), v, mask, Val{false}())
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


using VectorizationBase.LinearAlgebra: stride1
for v ∈ (:Vec, :SVec, :Val)
    vargs = Union{Symbol,Expr}[v === :Val ? :(::Val{W}) : :(::Type{$v{W,T}}), :(ptr::VectorizationBase.SparseStridedPointer{T})]
    # vcall = :(vmul(stride1(ptr), vrange(Val{W}())))
    vcall = :(vrangemul(Val{W}(), stride1(ptr), Val{0}()))
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
    scatter!(gep(gep(ptr, i), vrangemul(Val{W}(), stride1(ptr), Val{0}())), extract_data(v), Val{false}())
end
@inline function vstore!(ptr::VectorizationBase.AbstractSparseStridedPointer{T}, v::AbstractSIMDVector{W,T}, i, U::Unsigned) where {W,T}
    scatter!(gep(gep(ptr, i), vrangemul(Val{W}(), stride1(ptr), Val{0}())), extract_data(v), U, Val{false}())
end
@inline function vstore!(ptr::VectorizationBase.AbstractSparseStridedPointer{T}, v::AbstractSIMDVector{W,T}) where {W,T}
    scatter!(gep(ptr.ptr, vrangemul(Val{W}(), stride1(ptr), Val{0}())), extract_data(v), Val{false}())
end
@inline function vstore!(ptr::VectorizationBase.AbstractSparseStridedPointer{T}, v::AbstractSIMDVector{W,T}, U::Unsigned) where {W,T}
    scatter!(gep(ptr.ptr, vrangemul(Val{W}(), stride1(ptr), Val{0}())), extract_data(v), U, Val{false}())
end


# _MM support
# zero initialized
# scalar if only and Int
@inline vload(::AbstractZeroInitializedPointer{T}, ::Tuple{Int}) where {W,T<:Number} = zero(T)
# if a Vec of some kind, broadcast the zero
@inline vload(::AbstractZeroInitializedPointer{T}, ::Tuple{_MM{W},Vararg}) where {W,T<:Number} = vzero(SVec{W,T})
@inline vload(::AbstractZeroInitializedPointer{T}, ::Tuple{V,Vararg}) where {W,T<:Number,I<:Integer,V<:AbstractSIMDVector{W,I}} = vzero(SVec{W,T})
@inline vload(::AbstractZeroInitializedPointer{T}, ::Tuple{<:Integer,V,Vararg}) where {W,T<:Number,I<:Integer,V<:AbstractSIMDVector{W,I}} = vzero(SVec{W,T})
@inline vload(::AbstractZeroInitializedPointer{T}, ::Tuple{<:Integer,_MM{W},Vararg}) where {W,T<:Number} = vzero(SVec{W,T})
# peel off indices
@inline vload(ptr::AbstractZeroInitializedPointer{T}, i::Tuple{<:Integer,<:Integer,Vararg}) where {W,T<:Number} = vload(ptr, Base.tail(i))

vectypewidth(::Type{V}) where {W, V<:AbstractSIMDVector{W}} = W::Int
vectypewidth(::Type{_MM{W}}) where {W} = W::Int
vectypewidth(::Any) = 1
# returns expr for gep call, and bool-tuple (scalar,contiguous)
function packed_indexpr(Iparam, N, ::Type{T}) where {T}
    Ni = length(Iparam)
    Ni = min(Ni, N+1)
    Iₙ = Iparam[1]
    W::Int = vectypewidth(Iₙ)::Int
    indexpr = Expr(:ref, :i, 1)
    # if Iₙ <: _MM
    # end
    # check remaining indices.
    for n ∈ 2:Ni
        Iₙ = Iparam[n]
        Wₜ = vectypewidth(Iₙ)::Int
        W = W == 1 ? Wₜ : ((Wₜ == 1 || W == Wₜ) ? W : throw("$W ≠ $Wₜ but all vectors should be of the same width."))
        iexpr = Expr(:ref, :i, n)
        if Iₙ <: _MM
            # iexpr = Expr(:call, :+, Expr(:call, :svrange, Expr(:call, Expr(:curly, :Val, W)), T), Expr(:(.), iexpr, :i))
            iexpr = Expr(:call, :svrange, iexpr)
        end
        indexpr = Expr(:call, :muladd, iexpr, Expr(:ref, :s, n-1), indexpr)
    end
    W, Expr(:macrocall, Symbol("@inbounds"), LineNumberNode(@__LINE__, @__FILE__), indexpr)
end
function packed_strided_ptr_index(Iparam, N, ::Type{T}) where {T}
    W, indexpr = packed_indexpr(Iparam, N, T)
    W, Expr(:call, :gep, Expr(:(.), :ptr, QuoteNode(:ptr)), Expr(:call, :extract_data, indexpr))
end

@generated function vload(ptr::PackedStridedPointer{T,N}, i::I) where {T<:Number,N,I<:Tuple}
    W, gepcall = packed_strided_ptr_index(I.parameters, N, T)
    sexpr = Expr(:(=), :s, Expr(:(.), :ptr, QuoteNode(:strides)))
    if W == 1
        Expr(:block, Expr(:meta,:inline), sexpr, Expr(:call, :load, gepcall))
    else
        Expr(:block, Expr(:meta,:inline), sexpr, Expr(:call, :vload, Expr(:call, Expr(:curly, :Val, W)), gepcall))
    end
end
@generated function vload(ptr::PackedStridedPointer{T,N}, i::I, mask::Unsigned) where {T<:Number,N,I<:Tuple}
    W, gepcall = packed_strided_ptr_index(I.parameters, N, T)
    sexpr = Expr(:(=), :s, Expr(:(.), :ptr, QuoteNode(:strides)))
    if W == 1
        Expr(:block, Expr(:meta,:inline), sexpr, Expr(:call, :load, gepcall))
    else
        Expr(:block, Expr(:meta,:inline), sexpr, Expr(:call, :vload, Expr(:call, Expr(:curly, :Val, W)), gepcall, :mask))
    end
end
@generated function vstore!(ptr::AbstractPackedStridedPointer{T,N}, v::AbstractSIMDVector{1,T}, i::I) where {T<:Number,N,I<:Tuple}
    W, gepcall = packed_strided_ptr_index(I.parameters, N, T)
    sexpr = Expr(:(=), :s, Expr(:(.), :ptr, QuoteNode(:strides)))
    Expr(:block, Expr(:meta,:inline), sexpr, Expr(:call, :store!, gepcall, Expr(:macrocall, Symbol("@inbounds"), LineNumberNode(@__LINE__, @__FILE__), Expr(:call, :firstval, :v))))
end
@generated function vstore!(ptr::AbstractPackedStridedPointer{T,N}, v::AbstractSIMDVector{1,T}, i::I, mask::Unsigned) where {T<:Number,N,I<:Tuple}
    W, gepcall = packed_strided_ptr_index(I.parameters, N, T)
    sexpr = Expr(:(=), :s, Expr(:(.), :ptr, QuoteNode(:strides)))
    Expr(:block, Expr(:meta,:inline), sexpr, Expr(:call, :store!, gepcall, Expr(:macrocall, Symbol("@inbounds"), LineNumberNode(@__LINE__, @__FILE__), Expr(:call, :firstval, :v))))
end
@generated function vstore!(ptr::AbstractPackedStridedPointer{T,N}, v::AbstractSIMDVector{W1,T}, i::I) where {W1,T<:Number,N,I<:Tuple}
    W, gepcall = packed_strided_ptr_index(I.parameters, N, T)
    sexpr = Expr(:(=), :s, Expr(:(.), :ptr, QuoteNode(:strides)))
    Expr(:block, Expr(:meta,:inline), sexpr, Expr(:call, :vstore!, gepcall, Expr(:call, :extract_data, :v)))
end
@generated function vstore!(ptr::AbstractPackedStridedPointer{T,N}, v::AbstractSIMDVector{W1,T}, i::I, mask::Unsigned) where {W1,T<:Number,N,I<:Tuple}
    W, gepcall = packed_strided_ptr_index(I.parameters, N, T)
    sexpr = Expr(:(=), :s, Expr(:(.), :ptr, QuoteNode(:strides)))
    Expr(:block, Expr(:meta,:inline), sexpr, Expr(:call, :vstore!, gepcall, Expr(:call, :extract_data, :v), :mask))
end
@generated function vstore!(ptr::AbstractPackedStridedPointer{T,N}, sc::T, i::I, mask::Unsigned) where {T<:Number,N,I<:Tuple}
    W, gepcall = packed_strided_ptr_index(I.parameters, N, T)
    sexpr = Expr(:(=), :s, Expr(:(.), :ptr, QuoteNode(:strides)))
    if W == 1
        Expr(:block, Expr(:meta,:inline), sexpr, Expr(:call, :store!, gepcall, :sc))
    else
        Expr(:block, Expr(:meta,:inline), sexpr, Expr(:call, :vstore!, gepcall, Expr(:call, :vbroadcast, Expr(:curly, :Vec, W, T), :sc), :mask))
    end
end
@generated function vstore!(ptr::AbstractPackedStridedPointer{T,N}, sc::T, i::I) where {T<:Number,N,I<:Tuple}
    W, gepcall = packed_strided_ptr_index(I.parameters, N, T)
    sexpr = Expr(:(=), :s, Expr(:(.), :ptr, QuoteNode(:strides)))
    if W == 1
        Expr(:block, Expr(:meta,:inline), sexpr, Expr(:call, :store!, gepcall, :sc))
    else
        Expr(:block, Expr(:meta,:inline), sexpr, Expr(:call, :vstore!, gepcall, Expr(:call, :vbroadcast, Expr(:curly, :Vec, W, T), :sc)))
    end
end


using VectorizationBase: RowMajorStridedPointer, AbstractRowMajorStridedPointer
function rowmajor_strided_ptr_index(Iparam, N, ::Type{T}) where {T}
    Ni = length(Iparam)
    N = min(Ni - 1, N)
    Iₙ = Iparam[1]
    W::Int = vectypewidth(Iₙ)::Int
    indexpr = Expr(:ref, :i, N + 1)
    # check remaining indices.
    for n ∈ 1:N
        Iₙ = Iparam[n]
        Wₜ = vectypewidth(Iₙ)::Int
        W = W == 1 ? Wₜ : ((Wₜ == 1 || W == Wₜ) ? W : throw("$W ≠ $Wₜ but all vectors should be of the same width."))
        iexpr = Expr(:ref, :i, N + 1 - n)
        if Iₙ <: _MM
            # iexpr = Expr(:call, :+, Expr(:call, :svrange, Expr(:call, Expr(:curly, :Val, W)), T), Expr(:(.), iexpr, :i))
            iexpr = Expr(:call, :svrange, iexpr)
        end
        indexpr = Expr(:call, :muladd, iexpr, Expr(:ref, :s, n), indexpr)
    end
    indexpr = Expr(:macrocall, Symbol("@inbounds"), LineNumberNode(@__LINE__, @__FILE__), indexpr)
    W, Expr(:call, :gep, Expr(:(.), :ptr, QuoteNode(:ptr)), Expr(:call, :extract_data, indexpr))
end


@generated function vload(ptr::RowMajorStridedPointer{T,N}, i::I) where {T<:Number,N,I<:Tuple}
    W, gepcall = rowmajor_strided_ptr_index(I.parameters, N, T)
    sexpr = Expr(:(=), :s, Expr(:(.), :ptr, QuoteNode(:strides)))
    if W == 1
        Expr(:block, Expr(:meta,:inline), sexpr, Expr(:call, :load, gepcall))
    else
        Expr(:block, Expr(:meta,:inline), sexpr, Expr(:call, :vload, Expr(:call, Expr(:curly, :Val, W)), gepcall))
    end
end
@generated function vload(ptr::RowMajorStridedPointer{T,N}, i::I, mask::Unsigned) where {T<:Number,N,I<:Tuple}
    W, gepcall = rowmajor_strided_ptr_index(I.parameters, N, T)
    sexpr = Expr(:(=), :s, Expr(:(.), :ptr, QuoteNode(:strides)))
    if W == 1
        Expr(:block, Expr(:meta,:inline), sexpr, Expr(:call, :load, gepcall))
    else
        Expr(:block, Expr(:meta,:inline), sexpr, Expr(:call, :vload, Expr(:call, Expr(:curly, :Val, W)), gepcall, :mask))
    end
end
@generated function vstore!(ptr::AbstractRowMajorStridedPointer{T,N}, v::AbstractSIMDVector{1,T}, i::I) where {T<:Number,N,I<:Tuple}
    W, gepcall = rowmajor_strided_ptr_index(I.parameters, N, T)
    sexpr = Expr(:(=), :s, Expr(:(.), :ptr, QuoteNode(:strides)))
    Expr(:block, Expr(:meta,:inline), sexpr, Expr(:call, :store!, gepcall, Expr(:macrocall, Symbol("@inbounds"), LineNumberNode(@__LINE__, @__FILE__), Expr(:call, :firstval, :v))))
end
@generated function vstore!(ptr::AbstractRowMajorStridedPointer{T,N}, v::AbstractSIMDVector{1,T}, i::I, mask::Unsigned) where {T<:Number,N,I<:Tuple}
    W, gepcall = rowmajor_strided_ptr_index(I.parameters, N, T)
    sexpr = Expr(:(=), :s, Expr(:(.), :ptr, QuoteNode(:strides)))
    Expr(:block, Expr(:meta,:inline), sexpr, Expr(:call, :store!, gepcall, Expr(:macrocall, Symbol("@inbounds"), LineNumberNode(@__LINE__, @__FILE__), Expr(:call, :firstval, :v))))
end
@generated function vstore!(ptr::AbstractRowMajorStridedPointer{T,N}, v::AbstractSIMDVector{W1,T}, i::I) where {W1,T<:Number,N,I<:Tuple}
    W, gepcall = rowmajor_strided_ptr_index(I.parameters, N, T)
    sexpr = Expr(:(=), :s, Expr(:(.), :ptr, QuoteNode(:strides)))
    Expr(:block, Expr(:meta,:inline), sexpr, Expr(:call, :vstore!, gepcall, Expr(:call, :extract_data, :v)))
end
@generated function vstore!(ptr::AbstractRowMajorStridedPointer{T,N}, v::AbstractSIMDVector{W1,T}, i::I, mask::Unsigned) where {W1,T<:Number,N,I<:Tuple}
    W, gepcall = rowmajor_strided_ptr_index(I.parameters, N, T)
    sexpr = Expr(:(=), :s, Expr(:(.), :ptr, QuoteNode(:strides)))
    Expr(:block, Expr(:meta,:inline), sexpr, Expr(:call, :vstore!, gepcall, Expr(:call, :extract_data, :v), :mask))
end
@generated function vstore!(ptr::AbstractRowMajorStridedPointer{T,N}, sc::T, i::I, mask::Unsigned) where {T<:Number,N,I<:Tuple}
    W, gepcall = rowmajor_strided_ptr_index(I.parameters, N, T)
    sexpr = Expr(:(=), :s, Expr(:(.), :ptr, QuoteNode(:strides)))
    if W == 1
        Expr(:block, Expr(:meta,:inline), sexpr, Expr(:call, :store!, gepcall, :sc))
    else
        Expr(:block, Expr(:meta,:inline), sexpr, Expr(:call, :vstore!, gepcall, Expr(:call, :vbroadcast, Expr(:curly, :Vec, W, T), :sc), :mask))
    end
end
@generated function vstore!(ptr::AbstractRowMajorStridedPointer{T,N}, sc::T, i::I) where {T<:Number,N,I<:Tuple}
    W, gepcall = rowmajor_strided_ptr_index(I.parameters, N, T)
    sexpr = Expr(:(=), :s, Expr(:(.), :ptr, QuoteNode(:strides)))
    if W == 1
        Expr(:block, Expr(:meta,:inline), sexpr, Expr(:call, :store!, gepcall, :sc))
    else
        Expr(:block, Expr(:meta,:inline), sexpr, Expr(:call, :vstore!, gepcall, Expr(:call, :vbroadcast, Expr(:curly, :Vec, W, T), :sc)))
    end
end


using VectorizationBase: SparseStridedPointer, AbstractSparseStridedPointer
@inline VectorizationBase.gep(ptr::AbstractSparseStridedPointer, i::NTuple{W,Core.VecElement{I}}) where {W,I<:Integer} = gep(ptr.ptr, vmul(i,vbroadcast(Val{W}(), first(ptr.strides))))
function sparse_strided_ptr_index(Iparam, N, ::Type{T}) where {T}
    Ni = length(Iparam)
    Ni = N = min(Ni, N)
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
            # iexpr = Expr(:call, :+, Expr(:call, :svrange, Expr(:call, Expr(:curly, :Val, W)), T), iexpr)
            iexpr = Expr(:call, :svrange, iexpr)
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
@generated function vload(ptr::SparseStridedPointer{T,N}, i::I) where {T<:Number,N,I<:Tuple}
    W, gepcall = sparse_strided_ptr_index(I.parameters, N, T)
    sexpr = Expr(:(=), :s, Expr(:(.), :ptr, QuoteNode(:strides)))
    if W == 1
        Expr(:block, Expr(:meta,:inline), sexpr, Expr(:call, :load, gepcall))
    else
        Expr(:block, Expr(:meta,:inline), sexpr, Expr(:call, :vload, Expr(:call, Expr(:curly, :Val, W)), gepcall))
    end
end
@generated function vload(ptr::SparseStridedPointer{T,N}, i::I, mask::Unsigned) where {T<:Number,N,I<:Tuple}
    W, gepcall = sparse_strided_ptr_index(I.parameters, N, T)
    sexpr = Expr(:(=), :s, Expr(:(.), :ptr, QuoteNode(:strides)))
    if W == 1
        Expr(:block, Expr(:meta,:inline), sexpr, Expr(:call, :load, gepcall))
    else
        Expr(:block, Expr(:meta,:inline), sexpr, Expr(:call, :vload, Expr(:call, Expr(:curly, :Val, W)), gepcall, :mask))
    end
end
@generated function vstore!(ptr::AbstractSparseStridedPointer{T,N}, v::AbstractSIMDVector{1,T}, i::I) where {T<:Number,N,I<:Tuple}
    W, gepcall = sparse_strided_ptr_index(I.parameters, N, T)
    sexpr = Expr(:(=), :s, Expr(:(.), :ptr, QuoteNode(:strides)))
    Expr(:block, Expr(:meta,:inline), sexpr, Expr(:call, :store!, gepcall, Expr(:macrocall, Symbol("@inbounds"), LineNumberNode(@__LINE__, @__FILE__), Expr(:call, :firstval, :v))))
end
@generated function vstore!(ptr::AbstractSparseStridedPointer{T,N}, v::AbstractSIMDVector{1,T}, i::I, mask::Unsigned) where {T<:Number,N,I<:Tuple}
    W, gepcall = sparse_strided_ptr_index(I.parameters, N, T)
    sexpr = Expr(:(=), :s, Expr(:(.), :ptr, QuoteNode(:strides)))
    Expr(:block, Expr(:meta,:inline), sexpr, Expr(:call, :store!, gepcall, Expr(:macrocall, Symbol("@inbounds"), LineNumberNode(@__LINE__, @__FILE__), Expr(:call, :firstval, :v))))
end
@generated function vstore!(ptr::AbstractSparseStridedPointer{T,N}, v::AbstractSIMDVector{W1,T}, i::I) where {W1,T<:Number,N,I<:Tuple}
    W, gepcall = sparse_strided_ptr_index(I.parameters, N, T)
    sexpr = Expr(:(=), :s, Expr(:(.), :ptr, QuoteNode(:strides)))
    Expr(:block, Expr(:meta,:inline), sexpr, Expr(:call, :vstore!, gepcall, Expr(:call, :extract_data, :v)))
end
@generated function vstore!(ptr::AbstractSparseStridedPointer{T,N}, v::AbstractSIMDVector{W1,T}, i::I, mask::Unsigned) where {W1,T<:Number,N,I<:Tuple}
    W, gepcall = sparse_strided_ptr_index(I.parameters, N, T)
    sexpr = Expr(:(=), :s, Expr(:(.), :ptr, QuoteNode(:strides)))
    Expr(:block, Expr(:meta,:inline), sexpr, Expr(:call, :vstore!, gepcall, Expr(:call, :extract_data, :v), :mask))
end
@generated function vstore!(ptr::AbstractSparseStridedPointer{T,N}, sc::T, i::I) where {T<:Number,N,I<:Tuple}
    W, gepcall = sparse_strided_ptr_index(I.parameters, N, T)
    sexpr = Expr(:(=), :s, Expr(:(.), :ptr, QuoteNode(:strides)))
    if W == 1
        Expr(:block, Expr(:meta,:inline), sexpr, Expr(:call, :store!, gepcall, :sc))
    else
        Expr(:block, Expr(:meta,:inline), sexpr, Expr(:call, :vstore!, gepcall, Expr(:call, :vbroadcast, Expr(:curly, :Vec, W, T), :sc)))
    end
end
@generated function vstore!(ptr::AbstractSparseStridedPointer{T,N}, sc::T, i::I, mask::Unsigned) where {T<:Number,N,I<:Tuple}
    W, gepcall = sparse_strided_ptr_index(I.parameters, N, T)
    sexpr = Expr(:(=), :s, Expr(:(.), :ptr, QuoteNode(:strides)))
    if W == 1
        Expr(:block, Expr(:meta,:inline), sexpr, Expr(:call, :store!, gepcall, :sc))
    else
        Expr(:block, Expr(:meta,:inline), sexpr, Expr(:call, :vstore!, gepcall, Expr(:call, :vbroadcast, Expr(:curly, :Vec, W, T), :sc), :mask))
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
    N = min(length(Iparam), length(Xparam))
    Iₙ = Iparam[1]
    W::Int = vectypewidth(Iₙ)::Int
    # indexpr = Expr(:ref, :i, 1)
    local indexpr::Expr
    # check remaining indices.
    for n ∈ 1:N
        Iₙ = Iparam[n]
        Xₙ = (Xparam[n])::Int
        Wₜ = vectypewidth(Iₙ)::Int
        W = W == 1 ? Wₜ : ((Wₜ == 1 || W == Wₜ) ? W : throw("$W ≠ $Wₜ but all vectors should be of the same width."))
        iexpr = Expr(:ref, :i, n)
        if Xₙ > 1 && Iₙ <: _MM
            # iexpr = Expr(:call, :+, Expr(:call, :svrange, Expr(:call, Expr(:curly, :Val, W)), T), iexpr)
            iexpr = Expr(:call, :svrange, iexpr)
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
    W, Expr(:macrocall, Symbol("@inbounds"), LineNumberNode(@__LINE__, @__FILE__), indexpr)
end
function static_strided_ptr_gepcall(Iparam, Xparam, ::Type{T}) where {T}
    W, indexpr = static_strided_ptr_index(Iparam, Xparam, T)
    W, Expr(:call, :gep, Expr(:(.), :ptr, QuoteNode(:ptr)), Expr(:call, :extract_data, indexpr))
end
@generated function vload(ptr::StaticStridedPointer{T,X}, i::I) where {T<:Number,X,I<:Tuple}
    W, gepcall = static_strided_ptr_gepcall(I.parameters, X.parameters, T)
    if W == 1
        Expr(:block, Expr(:meta,:inline), Expr(:call, :load, gepcall))
    else
        Expr(:block, Expr(:meta,:inline), Expr(:call, :vload, Expr(:call, Expr(:curly, :Val, W)), gepcall))
    end
end
@generated function vload(ptr::StaticStridedPointer{T,X}, i::I, mask::Unsigned) where {T<:Number,X,I<:Tuple}
    W, gepcall = static_strided_ptr_gepcall(I.parameters, X.parameters, T)
    if W == 1
        Expr(:block, Expr(:meta,:inline), Expr(:call, :load, gepcall))
    else
        Expr(:block, Expr(:meta,:inline), Expr(:call, :vload, Expr(:call, Expr(:curly, :Val, W)), gepcall, :mask))
    end
end
@generated function vstore!(ptr::AbstractStaticStridedPointer{T,X}, v::AbstractSIMDVector{1,T}, i::I) where {T<:Number,X,I<:Tuple}
    W, gepcall = static_strided_ptr_gepcall(I.parameters, X.parameters, T)
    Expr(:block, Expr(:meta,:inline), Expr(:call, :store!, gepcall, Expr(:macrocall, Symbol("@inbounds"), LineNumberNode(@__LINE__, @__FILE__), Expr(:call, :firstval, :v))))
end
@generated function vstore!(ptr::AbstractStaticStridedPointer{T,X}, v::AbstractSIMDVector{W1,T}, i::I) where {W1,T<:Number,X,I<:Tuple}
    W, gepcall = static_strided_ptr_gepcall(I.parameters, X.parameters, T)
    Expr(:block, Expr(:meta,:inline), Expr(:call, :vstore!, gepcall, Expr(:call, :extract_data, :v)))
end
@generated function vstore!(ptr::AbstractStaticStridedPointer{T,X}, v::AbstractSIMDVector{1,T}, i::I, mask::Unsigned) where {T<:Number,X,I<:Tuple}
    W, gepcall = static_strided_ptr_gepcall(I.parameters, X.parameters, T)
    Expr(:block, Expr(:meta,:inline), Expr(:call, :store!, gepcall, Expr(:macrocall, Symbol("@inbounds"), LineNumberNode(@__LINE__, @__FILE__),  Expr(:call, :firstval, :v))))
end
@generated function vstore!(ptr::AbstractStaticStridedPointer{T,X}, v::AbstractSIMDVector{W1,T}, i::I, mask::Unsigned) where {W1,T<:Number,X,I<:Tuple}
    W, gepcall = static_strided_ptr_gepcall(I.parameters, X.parameters, T)
    Expr(:block, Expr(:meta,:inline), Expr(:call, :vstore!, gepcall, Expr(:call, :extract_data, :v), :mask))
end
@generated function vstore!(ptr::AbstractStaticStridedPointer{T,X}, sc::T, i::I) where {T<:Number,X,I<:Tuple}
    W, gepcall = static_strided_ptr_gepcall(I.parameters, X.parameters, T)
    if W == 1
        Expr(:block, Expr(:meta,:inline), Expr(:call, :store!, gepcall, :sc))
    else
        Expr(:block, Expr(:meta,:inline), Expr(:call, :vstore!, gepcall, Expr(:call, :vbroadcast, Expr(:curly, :Vec, W, T), :sc)))
    end
end
@generated function vstore!(ptr::AbstractStaticStridedPointer{T,X}, sc::T, i::I, mask::Unsigned) where {T<:Number,X,I<:Tuple}
    W, gepcall = static_strided_ptr_gepcall(I.parameters, X.parameters, T)
    if W == 1
        Expr(:block, Expr(:meta,:inline), Expr(:call, :store!, gepcall, :sc))
    else
        Expr(:block, Expr(:meta,:inline), Expr(:call, :vstore!, gepcall, Expr(:call, :vbroadcast, Expr(:curly, :Vec, W, T), :sc), :mask))
    end
end

@inline vstore!(ptr::AbstractStridedPointer{T}, s::Number, i::Tuple) where {T} = vstore!(ptr, convert(T, s), i)
@inline vstore!(ptr::AbstractStridedPointer{T}, s::Number, i::Tuple, mask::Unsigned) where {T} = vstore!(ptr, convert(T, s), i, mask)
@inline vstore!(ptr::AbstractStridedPointer{T1}, v::AbstractSIMDVector{W,T2}, i::Tuple) where {W,T1,T2} = vstore!(ptr, vconvert(Vec{W,T1}, v), i)
@inline vstore!(ptr::AbstractStridedPointer{T1}, v::AbstractSIMDVector{W,T2}, i::Tuple, mask::Unsigned) where {W,T1,T2} = vstore!(ptr, vconvert(Vec{W,T1}, v), i, mask)

using VectorizationBase: StaticStridedStruct

function staticstruct_vload_fromtup(t)
    Expr(
        :block,
        Expr(:meta,:inline),
        Expr(:(=), :s, Expr(:(.), :ptr, QuoteNode(:ptr))),
        Expr(:(=), :offset, Expr(:call, :+, Expr(:(.), :i, QuoteNode(:i)), Expr(:(.), :ptr, QuoteNode(:offset)))),
        Expr(:call, :SVec, Expr(:macrocall, Symbol("@inbounds"), LineNumberNode(@__LINE__, @__FILE__), t))
    )
end
@generated function vload(ptr::StaticStridedStruct, i::_MM{W}) where {W}
    t = Expr(:tuple)
    for w ∈ 1:W
        ref = Expr(:ref, :s, Expr(:call, :+, :offset, w))
        push!(t.args, Expr(:call, Expr(:(.), :Core, QuoteNode(:VecElement)), ref))
    end
    staticstruct_vload_fromtup(t)
end
@generated function vload(ptr::StaticStridedStruct, i::AbstractSIMDVector{W}) where {W}
    t = Expr(:tuple)
    for w ∈ 1:W
        ref = Expr(:ref, :s, Expr(:call, :+, :offset, Expr(:call, Expr(:(.), :VectorizationBase, QuoteNode(:extract_value), :i, w))))
        push!(t.args, Expr(:call, Expr(:(.), :Core, QuoteNode(:VecElement)), ref))
    end
    staticstruct_vload_fromtup(t)
end
@inline vload(ptr::StaticStridedStruct{T}, i::Union{_MM{W},<:AbstractSIMDVector{W}}, mask::Unsigned) where {W,T} = vifelse(mask, vload(ptr, i), vbroadcast(Val{W}(), zero(T)))

@generated function vload(ptr::StaticStridedStruct{T,X,S}, i::I) where {T,X,S,I<:Tuple}
    W, indexpr = static_strided_ptr_index(Iparam, Xparam, T)
    if W == 1
        Expr(:block, Expr(:meta, :inline), :(@inbounds ptr.ptr[1 + $indexpr]))
    else
        Expr(:block, Expr(:meta, :inline), Expr(:call, :vload, :ptr, indexpr))
    end
end
@generated function vload(ptr::StaticStridedStruct{T,X,S}, i::I, mask::Unsigned) where {T,X,S,I<:Tuple}
    W, indexpr = static_strided_ptr_index(Iparam, Xparam, T)
    if W == 1
        :(@inbounds ptr.ptr[1 + $indexpr])
    else
        Expr(:block, Expr(:meta, :inline), Expr(:call, :vload, :ptr, indexpr, :mask))
    end
end

@inline vload(r::AbstractRange{T}, i::Tuple{_MM{W}}) where {W,T} = SVec(vadd(vrangemul(Val{W}(), step(r), Val{0}()), @inbounds r[i[1].i + 1]))
#vmuladd(svrange(Val{W}(), T), step(r), @inbounds r[i[1].i + 1])
@inline vload(r::UnitRange{T}, i::Tuple{_MM{W}}) where {W,T} = @inbounds(_MM{W}(r[i[1].i + 1]))
# @inline vload(r::UnitRange{T}, i::Tuple{_MM{W}}) where {W,T} = svrangeincr(Val{W}(), @inbounds(r[i[1].i + 1]), Val{0}())
# Ignore masks
@inline vload(r::AbstractRange{T}, i::Tuple{_MM{W}}, ::Unsigned) where {W,T} = SVec(vadd(vrangemul(Val{W}(), step(r), Val{0}()), @inbounds r[i[1].i + 1]))
@inline vload(r::UnitRange{T}, i::Tuple{_MM{W}}, ::Unsigned) where {W,T} = @inbounds(_MM{W}(r[i[1].i + 1]))
# @inline vload(r::UnitRange{T}, i::Tuple{_MM{W}}, ::Unsigned) where {W,T} = svrangeincr(Val{W}(), @inbounds(r[i[1].i + 1]), Val{0}())



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
        push!(q.args, Expr(:macrocall, Symbol("@inbounds"), LineNumberNode(@__LINE__, @__FILE__), extractblock))
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
                    :call, :store!,
                    Expr(:call, :gep, :bptr, Expr(:tuple, n-1)),
                    Expr(
                        :call, :vsum,
                        Expr(:macrocall, Symbol("@inbounds"), LineNumberNode(@__LINE__, @__FILE__), Expr(:ref, :v, n))
                    )
                )
            )
        end
    end
    q
end




