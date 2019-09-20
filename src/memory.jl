

function valloc(::Type{T}, N::Int, sz::Int) where T
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
function valloc(f, ::Type{T}, N::Int, sz::Int) where T
    mem = valloc(T, N, sz)
    @inbounds for i in 1:sz
        mem[i] = f(i)
    end
    mem
end

@generated function sub2ind(dims::NTuple{N}, I::NTuple{N}) where N
    ex = :(I[$N] - 1)
    for i = (N - 1):-1:1
        ex = :(I[$i] - 1 + dims[$i] * $ex)
    end
    quote
        $(Expr(:meta, :inline))
        $ex + 1
    end
end

@generated function vload(::Type{Vec{N,T}}, ptr::Ptr{T},
                          ::Val{Aligned} = Val{false}()) where {N,T,Aligned}
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
    if align > 0
        push!(flags, "align $align")
    end
    push!(instrs, "%ptr = inttoptr $ptyp %0 to $vtyp*")
    push!(instrs, "%res = load $vtyp, $vtyp* %ptr" * join(flags, ", "))
    push!(instrs, "ret $vtyp %res")
    quote
        $(Expr(:meta, :inline))
        Base.llvmcall($((join(decls, "\n"), join(instrs, "\n"))),
            Vec{N,T}, Tuple{Ptr{T}}, ptr)
    end
end
@inline function vload(::Type{Vec{N,T}}, ptr::VectorizationBase.Pointer{T}, ::Val{Aligned} = Val{false}()) where {N,T,Aligned}
    vload(Vec{N,T}, ptr.ptr, Val{Aligned}())
end
@inline function vload(::Type{SVec{N,T}}, ptr::VectorizationBase.Pointer{T}, ::Val{Aligned} = Val{false}()) where {N,T,Aligned}
    SVec(vload(Vec{N,T}, ptr.ptr, Val{Aligned}()))
end
@inline function vload(::Type{SVec{N,T1}}, ptr::VectorizationBase.Pointer{T2}, ::Val{Aligned} = Val{false}()) where {N,T1,T2,Aligned}
    convert(SVec{N,T1}, SVec(vload(Vec{N,T2}, ptr.ptr, Val{Aligned}())))
end

@inline vloada(::Type{Vec{N,T}}, ptr::Ptr{T}) where {N,T} =
    vload(Vec{N,T}, ptr, Val{true}())


@inline function vload(::Type{Vec{N,T}},
                       ptr::Ptr{T},
                       i::Integer,
                       ::Val{Aligned} = Val{false}()) where {N,T,Aligned}
    #TODO @boundscheck 1 <= i <= length(arr) - (N-1) || throw(BoundsError())
    vload(Vec{N,T}, ptr + (i - 1)*sizeof(T), Val{Aligned}())
end
@inline function vload(::Type{Vec{N,T}},
                       arr::AbstractArray{T,D},
                       i::Integer,
                       ::Val{Aligned} = Val{false}()) where {N,T,Aligned,D}
    #TODO @boundscheck 1 <= i <= length(arr) - (N-1) || throw(BoundsError())
    vload(Vec{N,T}, pointer(arr, i), Val{Aligned}())
end
@inline function vload(::Type{Vec{N,T}},
                       arr::AbstractArray{T,D},
                       ind::NTuple{D,<:Integer},
                       ::Val{Aligned} = Val{false}()) where {N,T,Aligned,D}
    #TODO @boundscheck 1 <= i <= length(arr) - (N-1) || throw(BoundsError())
    i = sub2ind(size(arr), ind)
    vload(Vec{N,T}, pointer(arr, i), Val{Aligned}())
end
@inline function vloada(::Type{Vec{N,T}},
                        arr::AbstractArray{T,1},
                        i::Integer) where {N,T}
    vload(Vec{N,T}, arr, i, Val{true}())
end



@inline vload(::Type{SVec{N,T}}, ptr::Ptr{T}, ::Val{Aligned} = Val{false}()) where {N,T,Aligned} =
    SVec(vload(Vec{N,T}, ptr, Val{Aligned}()))

@inline vloada(::Type{SVec{N,T}}, ptr::Ptr{T}) where {N,T} =
    SVec(vload(Vec{N,T}, ptr, Val{true}()))


@inline function vload(::Type{SVec{N,T}},
                       ptr::Ptr{T},
                       i::Integer,
                       ::Val{Aligned} = Val{false}()) where {N,T,Aligned}
    #TODO @boundscheck 1 <= i <= length(arr) - (N-1) || throw(BoundsError())
    SVec(vload(Vec{N,T}, ptr + (i - 1)*sizeof(T), Val{Aligned}()))
end
@inline function vload(::Type{SVec{N,T}},
                       arr::AbstractArray{T,D},
                       i::Integer,
                       ::Val{Aligned} = Val{false}()) where {N,T,Aligned,D}
    #TODO @boundscheck 1 <= i <= length(arr) - (N-1) || throw(BoundsError())
    SVec(vload(Vec{N,T}, pointer(arr, i), Val{Aligned}()))
end
@inline function vload(::Type{SVec{N,T}},
                       arr::AbstractArray{T,D},
                       ind::NTuple{D,<:Integer},
                       ::Val{Aligned} = Val{false}()) where {N,T,Aligned,D}
    #TODO @boundscheck 1 <= i <= length(arr) - (N-1) || throw(BoundsError())
    i = sub2ind(size(arr), ind)
    SVec(vload(Vec{N,T}, pointer(arr, i), Val{Aligned}()))
end
@inline function vloada(::Type{SVec{N,T}},
                        arr::AbstractArray{T,1},
                        i::Integer) where {N,T}
    SVec(vload(Vec{N,T}, arr, i, Val{true}()))
end


@generated function vload(::Type{Vec{N,T}}, ptr::Ptr{T},
                          mask::Vec{N,Bool},
                          ::Val{Aligned} = Val{false}()) where {N,T,Aligned}
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
            Vec{N,T}, Tuple{Ptr{T}, Vec{N,Bool}}, ptr, mask)
    end
end
@generated function vload(
    ::Type{Vec{N,T}}, ptr::Ptr{T}, mask::U, ::Val{Aligned} = Val{false}()
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
@inline function vload(::Type{Vec{N,T}}, ptr::VectorizationBase.Pointer{T},
                mask::Union{Vec{N,Bool},SVec{N,Bool},<:Unsigned}, ::Val{Aligned} = Val{false}()) where {N,T,Aligned}
    vload(Vec{N,T}, ptr.ptr, mask, Val{Aligned}())
end
@inline function vload(::Type{SVec{N,T}}, ptr::VectorizationBase.Pointer{T},
                mask::Union{Vec{N,Bool},SVec{N,Bool},<:Unsigned}, ::Val{Aligned} = Val{false}()) where {N,T,Aligned}
    SVec(vload(Vec{N,T}, ptr.ptr, mask, Val{Aligned}()))
end
@inline function vload(::Type{SVec{N,T1}}, ptr::VectorizationBase.Pointer{T2},
                mask::Union{Vec{N,Bool},SVec{N,Bool},<:Unsigned}, ::Val{Aligned} = Val{false}()) where {N,T1,T2,Aligned}
    convert(SVec{N,T1}, SVec(vload(Vec{N,T2}, ptr.ptr, mask, Val{Aligned}())))
end

@inline vloada(::Type{Vec{N,T}}, ptr::Ptr{T}, mask::Vec{N,Bool}) where {N,T} =
    vload(Vec{N,T}, ptr, mask, Val{true}())

@inline function vload(::Type{Vec{N,T}},
                       ptr::Ptr{T},
                       i::Integer, mask::Union{Vec{N,Bool},<:Unsigned},
                       ::Val{Aligned} = Val{false}()) where {N,T,Aligned}
    #TODO @boundscheck 1 <= i <= length(arr) - (N-1) || throw(BoundsError())
    vload(Vec{N,T}, ptr + (i - 1)*sizeof(T), mask, Val{Aligned}())
end
@inline function vload(::Type{Vec{N,T}},
                       arr::AbstractArray{T},
                       i::Integer, mask::Union{Vec{N,Bool},<:Unsigned},
                       ::Val{Aligned} = Val{false}()) where {N,T,Aligned}
    #TODO @boundscheck 1 <= i <= length(arr) - (N-1) || throw(BoundsError())
    vload(Vec{N,T}, pointer(arr, i), mask, Val{Aligned}())
end
@inline function vload(::Type{Vec{N,T}},
                       arr::AbstractArray{T,D},
                       ind::NTuple{D,<:Integer}, mask::Union{Vec{N,Bool},<:Unsigned},
                       ::Val{Aligned} = Val{false}()) where {N,T,D,Aligned}
    #TODO @boundscheck 1 <= i <= length(arr) - (N-1) || throw(BoundsError())
    i = sub2ind(size(arr), ind)
    vload(Vec{N,T}, pointer(arr, i), mask, Val{Aligned}())
end
@inline function vloada(::Type{Vec{N,T}},
                        arr::AbstractArray{T}, i::Integer,
                        mask::Union{Vec{N,Bool},<:Unsigned}) where {N,T}
    vload(Vec{N,T}, arr, i, mask, Val{true}())
end

@inline vload(::Type{SVec{N,T}}, ptr::Ptr{T}, mask::Union{Vec{N,Bool},<:Unsigned}, ::Val{Aligned} = Val{false}()) where {N,T,Aligned} =
    SVec(vload(Vec{N,T}, ptr, mask, Val{Aligned}()))

@inline vloada(::Type{SVec{N,T}}, ptr::Ptr{T}, mask::Union{Vec{N,Bool},<:Unsigned}) where {N,T} =
    SVec(vload(Vec{N,T}, ptr, mask, Val{true}()))

@inline function vload(::Type{SVec{N,T}},
                       ptr::Ptr{T},
                       i::Integer, mask::Union{Vec{N,Bool},<:Unsigned},
                       ::Val{Aligned} = Val{false}()) where {N,T,Aligned}
    #TODO @boundscheck 1 <= i <= length(arr) - (N-1) || throw(BoundsError())
    SVec(vload(Vec{N,T}, ptr + (i - 1)*sizeof(T), mask, Val{Aligned}()))
end
@inline function vload(::Type{SVec{N,T}},
                       arr::AbstractArray{T},
                       i::Integer, mask::Union{Vec{N,Bool},<:Unsigned},
                       ::Val{Aligned} = Val{false}()) where {N,T,Aligned}
    #TODO @boundscheck 1 <= i <= length(arr) - (N-1) || throw(BoundsError())
    SVec(vload(Vec{N,T}, pointer(arr, i), mask, Val{Aligned}()))
end
@inline function vload(::Type{SVec{N,T}},
                       arr::AbstractArray{T,D},
                       ind::NTuple{D,<:Integer}, mask::Union{Vec{N,Bool},<:Unsigned},
                       ::Val{Aligned} = Val{false}()) where {N,T,D,Aligned}
    #TODO @boundscheck 1 <= i <= length(arr) - (N-1) || throw(BoundsError())
    i = sub2ind(size(arr), ind)
    SVec(vload(Vec{N,T}, pointer(arr, i), mask, Val{Aligned}()))
end
@inline function vloada(::Type{SVec{N,T}},
                        arr::AbstractArray{T}, i::Integer,
                        mask::Union{Vec{N,Bool},<:Unsigned}) where {N,T}
    SVec(vload(Vec{N,T}, arr, i, mask, Val{true}()))
end


@generated function vstore!(ptr::Ptr{T}, v::Vec{N,T},
                           ::Val{Aligned} = Val{false}()) where {N,T,Aligned}
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
    if align > 0
        push!(flags, "align $align")
    end
    push!(instrs, "%ptr = inttoptr $ptyp %0 to $vtyp*")
    push!(instrs, "store $vtyp %1, $vtyp* %ptr" * join(flags, ", "))
    push!(instrs, "ret void")
    quote
        $(Expr(:meta, :inline))
        Base.llvmcall($((join(decls, "\n"), join(instrs, "\n"))),
                      Cvoid, Tuple{Ptr{T}, Vec{N,T}}, ptr, v)
    end
end
@inline function vstore!(ptr::Ptr{T}, v::AbstractStructVec{N,T}, ::Val{Aligned} = Val{false}()) where {N,T,Aligned}
    vstore!(ptr, extract_data(v), Val{Aligned}())
end

@inline vstorea!(ptr::Ptr{T}, v::AbstractSIMDVector{N,T}) where {N,T} =
            vstore!(ptr, extract_data(v), Val{true}())


@inline function vstore!(ptr::Ptr{T}, v::AbstractSIMDVector{N,T},
                        i::Integer,
                        ::Val{Aligned} = Val{false}()) where {N,T,Aligned}
    vstore!(ptr + (i - 1)*sizeof(T), extract_data(v), Val{Aligned}())
end
@inline function vstore!(arr::AbstractArray{T,D},
                        v::AbstractSIMDVector{N,T},
                        i::Integer,
                        ::Val{Aligned} = Val{false}()) where {N,T,Aligned,D}
    @boundscheck 1 <= i <= length(arr) - (N-1) || throw(BoundsError())
    vstore!(pointer(arr, i), extract_data(v), Val{Aligned}())
end
@inline function vstore!(arr::AbstractArray{T,D},
                        v::AbstractSIMDVector{N,T},
                        ind::NTuple,
                        ::Val{Aligned} = Val{false}()) where {N,T,Aligned,D}
    i = sub2ind(size(arr), ind)
    # @boundscheck 1 <= i <= length(arr) - (N-1) || throw(BoundsError())
    vstore!(pointer(arr, i), extract_data(v), Val{Aligned}())
end
@inline function vstorea!(arr::AbstractArray{T,D}, v::AbstractSIMDVector{N,T}, i::Integer) where {N,T,D}
    vstore!(arr, extract_data(v), i, Val{true}())
end

@generated function vstore!(ptr::Ptr{T}, v::Vec{N,T}, mask::Vec{N,Bool},
                           ::Val{Aligned} = Val{false}()) where {N,T,Aligned}
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
        Base.llvmcall($((join(decls, "\n"), join(instrs, "\n"))),
            Cvoid, Tuple{Ptr{T}, Vec{N,T}, Vec{N,Bool}},
            ptr, v, mask)
    end
end
@generated function vstore!(ptr::Ptr{T}, v::Vec{N,T}, mask::U,
                           ::Val{Aligned} = Val{false}()) where {N,T,Aligned,U<:Unsigned}
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
@inline function vstore!(ptr::Ptr{T}, v::AbstractStructVec{N,T}, mask::Union{Vec{N,Bool},Unsigned}, ::Val{Aligned} = Val{false}()) where {N,T,Aligned}
    vstore!(ptr, extract_data(v), mask, Val{Aligned}())
end

@inline vstorea!(ptr::Ptr{T}, v::AbstractSIMDVector{N,T}, mask::Union{Vec{N,Bool},Unsigned}) where {N,T} =
    vstore!(ptr, extract_data(v), mask, Val{true}())

@inline function vstore!(ptr::Ptr{T},
                        v::AbstractSIMDVector{N,T},
                        i::Integer,
                        mask::Union{Vec{N,Bool},Unsigned},
                        ::Val{Aligned} = Val{false}()) where {N,T,Aligned}
    #TODO @boundscheck 1 <= i <= length(arr) - (N-1) || throw(BoundsError())
    vstore!(ptr + (i - 1)*sizeof(T), extract_data(v), mask, Val{Aligned}())
end
@inline function vstore!(arr::AbstractArray{T,D},
                        v::AbstractSIMDVector{N,T},
                        i::Integer,
                        mask::Union{Vec{N,Bool},Unsigned},
                        ::Val{Aligned} = Val{false}()) where {N,T,Aligned,D}
    #TODO @boundscheck 1 <= i <= length(arr) - (N-1) || throw(BoundsError())
    vstore!(pointer(arr, i), extract_data(v), mask, Val{Aligned}())
end
@inline function vstorea!(arr::AbstractArray{T},
                         v::AbstractSIMDVector{N,T},
                         i::Integer, mask::Union{Vec{N,Bool},Unsigned}) where {N,T}
    vstore!(arr, extract_data(v), i, mask, Val{true}())
end


@inline function vstore!(ptr::VectorizationBase.Pointer{T}, v::Vec{N,T}, ::Val{Aligned} = Val{false}()) where {N,T,Aligned}
    vstore!(ptr.ptr, v, Val{Aligned}())
end
@inline function vstore!(ptr::VectorizationBase.Pointer{T}, v::AbstractSIMDVector{N,T}, ::Val{Aligned} = Val{false}()) where {N,T,Aligned}
    vstore!(ptr.ptr, extract_data(v), Val{Aligned}())
end
@inline function vstore!(ptr::VectorizationBase.Pointer{T}, v::Vec{N,T},
                mask::Union{Vec{N,Bool},SVec{N,Bool},<:Unsigned}, ::Val{Aligned} = Val{false}()) where {N,T,Aligned}
    vstore!(ptr.ptr, v, mask, Val{Aligned}())
end
@inline function vstore!(ptr::VectorizationBase.Pointer{T}, v::AbstractSIMDVector{N,T},
                mask::Union{Vec{N,Bool},SVec{N,Bool},<:Unsigned}, ::Val{Aligned} = Val{false}()) where {N,T,Aligned}
    vstore!(ptr.ptr, extract_data(v), mask, Val{Aligned}())
end



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
        "%res = call $vtyp @llvm.masked.gather.$(suffix(N,T)).$(suffix(N,Ptr{T}))($vptrtyp %ptr, i32 $align, <$N x i1> <i1 true$(mask)>, $vtyp $(llvmconst(N, T, 0)))")
    push!(instrs, "ret $vtyp %res")
    quote
        $(Expr(:meta, :inline))
        Base.llvmcall($((join(decls, "\n"), join(instrs, "\n"))),
            Vec{$N,$T}, Tuple{Vec{$N,Ptr{$T}}}, ptr)
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
        Base.llvmcall($((join(decls, "\n"), join(instrs, "\n"))),
            Vec{$N,$T}, Tuple{Vec{$N,Ptr{$T}}, $U}, ptr, mask)
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
    # push!(decls,
        # "declare void @llvm.masked.scatter.$(suffix(N,T)).$(suffix(N,Ptr{T}))($vtyp, $vptrtyp, i32, <$N x i1>)")
    # push!(instrs,
        # "call void @llvm.masked.scatter.$(suffix(N,T)).$(suffix(N,Ptr{T}))($vtyp %0, $vptrtyp %ptr, i32 $align, <$N x i1> <i1 true$(mask)>)")
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
        push!(instrs, "%masktrunc = trunc $mtyp_input %3 to $mtyp_trunc")
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
        Base.llvmcall($((join(decls, "\n"), join(instrs, "\n"))),
            Cvoid, Tuple{Vec{$N,$T}, Vec{$N,Ptr{$T}}, $U}, v, ptr, mask)
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
