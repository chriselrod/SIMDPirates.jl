@generated function vrangeincr(::Val{W}, i::I, ::Val{O}) where {W,I<:Integer,O}
    bytes = I === Int ? min(8, VectorizationBase.prevpow2(VectorizationBase.REGISTER_SIZE ÷ W)) : sizeof(I)
    # bytes = min(8, VectorizationBase.prevpow2(VectorizationBase.REGISTER_SIZE ÷ W))
    bits = 8bytes
    jtypesym = Symbol(:Int, bits)
    iexpr = bytes == sizeof(I) ? :i : Expr(:call, :%, :i, jtypesym)
    typ = "i$(bits)"
    vtyp = "<$W x $typ>"
    rangevec = join(("$typ $(w+O)" for w ∈ 0:W-1), ", ")
    instrs = String[]
    push!(instrs, "%ie = insertelement $vtyp undef, $typ %0, i32 0")
    push!(instrs, "%v = shufflevector $vtyp %ie, $vtyp undef, <$W x i32> zeroinitializer")
    push!(instrs, "%res = add nsw $vtyp %v, <$rangevec>")
    push!(instrs, "ret $vtyp %res")
    quote
        $(Expr(:meta,:inline))
        Base.llvmcall(
            $(join(instrs,"\n")), Vec{$W,$jtypesym}, Tuple{$jtypesym}, $iexpr
        )
    end
end
@generated function vrangeincr(::Val{W}, i::T, ::Val{O}) where {W,T<:FloatingTypes,O}
    typ = llvmtype(T)
    vtyp = "<$W x $typ>"
    rangevec = join(("$typ $(w+O).0" for w ∈ 0:W-1), ", ")
    instrs = String[]
    push!(instrs, "%ie = insertelement $vtyp undef, $typ %0, i32 0")
    push!(instrs, "%v = shufflevector $vtyp %ie, $vtyp undef, <$W x i32> zeroinitializer")
    push!(instrs, "%res = fadd $vtyp %v, <$rangevec>")
    push!(instrs, "ret $vtyp %res")
    quote
        $(Expr(:meta,:inline))
        Base.llvmcall(
            $(join(instrs,"\n")), Vec{$W,$T}, Tuple{$T}, i
        )
    end
end
@generated function vrangemul(::Val{W}, i::I, ::Val{O}) where {W,I<:Integer,O}
    bytes = I === Int ? min(8, VectorizationBase.prevpow2(VectorizationBase.REGISTER_SIZE ÷ W)) : sizeof(I)
    bits = 8bytes
    jtypesym = Symbol(:Int, bits)
    iexpr = bytes == sizeof(I) ? :i : Expr(:call, :%, :i, jtypesym)
    typ = "i$(bits)"
    vtyp = "<$W x $typ>"
    rangevec = join(("$typ $(w+O)" for w ∈ 0:W-1), ", ")
    instrs = String[]
    push!(instrs, "%ie = insertelement $vtyp undef, $typ %0, i32 0")
    push!(instrs, "%v = shufflevector $vtyp %ie, $vtyp undef, <$W x i32> zeroinitializer")
    push!(instrs, "%res = mul nsw $vtyp %v, <$rangevec>")
    push!(instrs, "ret $vtyp %res")
    quote
        $(Expr(:meta,:inline))
        Base.llvmcall(
            $(join(instrs,"\n")), Vec{$W,$jtypesym}, Tuple{$jtypesym}, $iexpr
        )
    end
end
@generated function vrangemul(::Val{W}, i::T, ::Val{O}) where {W,T<:FloatingTypes,O}
    typ = llvmtype(T)
    vtyp = "<$W x $typ>"
    rangevec = join(("$typ $(w+O).0" for w ∈ 0:W-1), ", ")
    instrs = String[]
    push!(instrs, "%ie = insertelement $vtyp undef, $typ %0, i32 0")
    push!(instrs, "%v = shufflevector $vtyp %ie, $vtyp undef, <$W x i32> zeroinitializer")
    push!(instrs, "%res = fmul fast $vtyp %v, <$rangevec>")
    push!(instrs, "ret $vtyp %res")
    quote
        $(Expr(:meta,:inline))
        Base.llvmcall(
            $(join(instrs,"\n")), Vec{$W,$T}, Tuple{$T}, i
        )
    end
end
@inline svrangeincr(::Val{W}, i, ::Val{O}) where {W,O} = SVec(vrangeincr(Val{W}(), i, Val{O}()))
@inline svrangemul(::Val{W}, i, ::Val{O}) where {W,O} = SVec(vrangemul(Val{W}(), i, Val{O}()))


@inline vrange(i::_MM{W}) where {W} = vrangeincr(Val{W}(), i.i, Val{0}())
@inline svrange(i::_MM{W}) where {W} = SVec(vrangeincr(Val{W}(), i.i, Val{0}()))
@inline Base.:(+)(i::_MM{W}, j::_MM{W}) where {W} = SVec(vadd(vrange(i), vrange(j)))
@inline Base.:(+)(i::_MM{W}, j::AbstractSIMDVector{W}) where {W} = vadd(vrange(i), j)
@inline Base.:(+)(i::AbstractSIMDVector{W}, j::_MM{W}) where {W} = vadd(i, vrange(j))
@inline Base.:(*)(i::_MM{W}, j::AbstractSIMDVector{W}) where {W} = vmul(vrange(i), j)
@inline Base.:(*)(i::AbstractSIMDVector{W}, j::_MM{W}) where {W} = vmul(i, vrange(j))
@inline vadd(i::_MM{W}, j::_MM{W}) where {W} = SVec(vadd(vrange(i), vrange(j)))
@inline vadd(i::_MM{W}, j::AbstractSIMDVector{W}) where {W} = vadd(vrange(i), j)
@inline vadd(i::AbstractSIMDVector{W}, j::_MM{W}) where {W} = vadd(i, vrange(j))
@inline vmul(i::_MM{W}, j::AbstractSIMDVector{W}) where {W} = vmul(vrange(i), j)
@inline vmul(i::AbstractSIMDVector{W}, j::_MM{W}) where {W} = vmul(i, vrange(j))


@inline vrange(::Val{W}) where {W} = vrange(Val{W}(), Float64)
@inline svrange(::Val{W}) where {W} = svrange(Val{W}(), Float64)

@inline vrange(i::_MM{W}, ::Type{T}) where {W,T} = vrangeincr(Val{W}(), T(i.i), Val{0}())
@inline vrange(i::_MM{W}, ::Type{T}) where {W,T <: Integer} = vrangeincr(Val{W}(), i.i % T, Val{0}())
@inline svrange(i::_MM, ::Type{T}) where {T} = SVec(vrange(i, T))


@inline Base.:(<<)(i::_MM, j::Integer) = svrange(i) << j
@inline Base.:(>>)(i::_MM, j::Integer) = svrange(i) >> j
@inline Base.:(>>>)(i::_MM, j::Integer) = svrange(i) >>> j

@inline Base.:(*)(i::_MM{W}, j::T) where {W,T} = vmul(svrange(i), j)
@inline Base.:(*)(j::T, i::_MM{W}) where {W,T} = vmul(svrange(i), j)
@inline vmul(i::_MM{W}, j::T) where {W,T} = vmul(svrange(i), j)
@inline vmul(j::T, i::_MM{W}) where {W,T} = vmul(svrange(i), j)
@inline vconvert(::Type{Vec{W,T}}, i::_MM{W}) where {W,T} = vrange(i, T)
@inline vconvert(::Type{SVec{W,T}}, i::_MM{W}) where {W,T} = svrange(i, T)




@inline Base.:(-)(i::Integer, j::_MM{W}) where {W} = i - svrange(j)
@inline Base.:(-)(::Static{i}, j::_MM{W}) where {W,i} = i - svrange(j)
@inline Base.:(-)(i::_MM{W}, j::_MM{W}) where {W} = svrange(i) - svrange(j)
@inline Base.:(-)(i::_MM{W}) where {W} = -svrange(i)
@inline vsub(i::Integer, j::_MM{W}) where {W} = i - svrange(j)
@inline vsub(::Static{i}, j::_MM{W}) where {W,i} = i - svrange(j)
@inline vsub(i::_MM{W}, j::_MM{W}) where {W} = svrange(i) - svrange(j)
@inline vsub(i::_MM{W}) where {W} = -svrange(i)



@inline Base.:(<)(i::_MM, j::Integer) = svrange(i) < j
@inline Base.:(<)(i::Integer, j::_MM) = i < svrange(j)
@inline Base.:(<)(i::_MM, ::Static{j}) where {j} = svrange(i) < j
@inline Base.:(<)(::Static{i}, j::_MM) where {i} = i < svrange(j)
@inline Base.:(<)(i::_MM, j::_MM) = svrange(i) < svrange(j)
@inline Base.:(>)(i::_MM, j::Integer) = svrange(i) > j
@inline Base.:(>)(i::Integer, j::_MM) = i > svrange(j)
@inline Base.:(>)(i::_MM, ::Static{j}) where {j} = svrange(i) > j
@inline Base.:(>)(::Static{i}, j::_MM) where {i} = i > svrange(j)
@inline Base.:(>)(i::_MM, j::_MM) = svrange(i) > svrange(j)
@inline Base.:(==)(i::_MM, j::Integer) = svrange(i) == j
@inline Base.:(==)(i::Integer, j::_MM) = i == svrange(j)
@inline Base.:(==)(i::_MM, ::Static{j}) where {j} = svrange(i) == j
@inline Base.:(==)(::Static{i}, j::_MM) where {i} = i == svrange(j)
@inline Base.:(==)(i::_MM, j::_MM) = svrange(i) == svrange(j)
@inline Base.:(!=)(i::_MM, j::Integer) = svrange(i) != j
@inline Base.:(!=)(i::Integer, j::_MM) = i != svrange(j)
@inline Base.:(!=)(i::_MM, ::Static{j}) where {j} = svrange(i) != j
@inline Base.:(!=)(::Static{i}, j::_MM) where {i} = i != svrange(j)
@inline Base.:(!=)(i::_MM, j::_MM) = svrange(i) != svrange(j)
@inline Base.:(&)(i::_MM, j::Integer) = svrange(i) & j
@inline Base.:(&)(i::Integer, j::_MM) = i & svrange(j)
@inline Base.:(&)(i::_MM, ::Static{j}) where {j} = svrange(i) & j
@inline Base.:(&)(::Static{i}, j::_MM) where {i} = i & svrange(j)
@inline Base.:(&)(i::_MM, j::_MM) = svrange(i) & svrange(j)
@inline Base.:(|)(i::_MM, j::Integer) = svrange(i) | j
@inline Base.:(|)(i::Integer, j::_MM) = i | svrange(j)
@inline Base.:(|)(i::_MM, ::Static{j}) where {j} = svrange(i) | j
@inline Base.:(|)(::Static{i}, j::_MM) where {i} = i | svrange(j)
@inline Base.:(|)(i::_MM, j::_MM) = svrange(i) | svrange(j)
@inline Base.:(⊻)(i::_MM, j::Integer) = svrange(i) ⊻ j
@inline Base.:(⊻)(i::Integer, j::_MM) = i ⊻ svrange(j)
@inline Base.:(⊻)(i::_MM, ::Static{j}) where {j} = svrange(i) ⊻ j
@inline Base.:(⊻)(::Static{i}, j::_MM) where {i} = i ⊻ svrange(j)
@inline Base.:(⊻)(i::_MM, j::_MM) = svrange(i) ⊻ svrange(j)
@inline Base.:(*)(i::_MM, j::_MM) = SVec(vmul(vrange(i), vrange(j)))
@inline vmul(i::_MM, j::_MM) = SVec(vmul(vrange(i), vrange(j)))


