function tuple_range_vector_expr(W)
    t = Expr(:tuple)
    for w ∈ zero(W):W-one(W)
        push!(t.args, Expr(:call, Expr(:(.), :Core, QuoteNode(:VecElement)), w))
    end
    t
end
@generated function vrangeincr(::Val{W}, i::I, ::Val{O}) where {W,I<:Integer,O}
    bytes = min(8, VectorizationBase.prevpow2(VectorizationBase.REGISTER_SIZE ÷ W))
    bits = 8bytes
    jtypesym = Symbol(:Int, bits)
    iexpr = bytes == sizeof(I) ? :i : Expr(:call, :%, :i, jtypesym)
    typ = "i$(bits)"
    vtyp = "<$W x $typ>"
    rangevec = join(("$typ $(w+O)" for w ∈ 0:W-1), ", ")
    instrs = String[]
    push!(instrs, "%ie = insertelement $vtyp undef, $typ %0, i32 0")
    push!(instrs, "%v = shufflevector $vtyp %ie, $vtyp undef, <$W x i32> zeroinitializer")
    push!(instrs, "%res = add $vtyp %v, <$rangevec>")
    push!(instrs, "ret $vtyp %res")
    quote
        $(Expr(:meta,:inline))
        Base.llvmcall(
            $(join(instrs,"\n")), Vec{$W,$jtypesym}, Tuple{$jtypesym}, $iexpr
        )
    end
end
@generated function vrangemul(::Val{W}, i::I, ::Val{O}) where {W,I<:Integer,O}
    bytes = min(8, VectorizationBase.prevpow2(VectorizationBase.REGISTER_SIZE ÷ W))
    bits = 8bytes
    jtypesym = Symbol(:Int, bits)
    iexpr = bytes == sizeof(I) ? :i : Expr(:call, :%, :i, jtypesym)
    typ = "i$(bits)"
    vtyp = "<$W x $typ>"
    rangevec = join(("$typ $(w+O)" for w ∈ 0:W-1), ", ")
    instrs = String[]
    push!(instrs, "%ie = insertelement $vtyp undef, $typ %0, i32 0")
    push!(instrs, "%v = shufflevector $vtyp %ie, $vtyp undef, <$W x i32> zeroinitializer")
    push!(instrs, "%res = mul $vtyp %v, <$rangevec>")
    push!(instrs, "ret $vtyp %res")
    quote
        $(Expr(:meta,:inline))
        Base.llvmcall(
            $(join(instrs,"\n")), Vec{$W,$jtypesym}, Tuple{$jtypesym}, $iexpr
        )
    end
end
@inline svrangeincr(::Val{W}, i, ::Val{O}) where {W,O} = SVec(vrangeincr(Val{W}(), i, Val{O}()))
@inline svrangemul(::Val{W}, i, ::Val{O}) where {W,O} = SVec(vrangemul(Val{W}(), i, Val{O}()))
function intrangetuple(W, ::Type{T}) where {T}
    ret::Expr = if sizeof(T) == 8
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
    ret
end
@generated function vrange(::Val{W}, ::Type{T}) where {W, T}
    Expr(:block, Expr(:meta,:inline), intrangetuple(W, T))
end
@generated function svrange(::Val{W}, ::Type{T}) where {W, T}
    Expr(:block, Expr(:meta,:inline), Expr(:call, :SVec, intrangetuple(W, T)))
end


using VectorizationBase: _MM, AbstractZeroInitializedPointer
@inline vrange(i::_MM{W}) where {W} = vrangeincr(Val{W}(), i.i, Val{0}())
@inline svrange(i::_MM{W}) where {W} = SVec(vrangeincr(Val{W}(), i.i, Val{0}()))
@inline Base.:(+)(i::_MM{W}, j::AbstractSIMDVector{W}) where {W} = vadd(vrange(i), j)
@inline Base.:(+)(i::AbstractSIMDVector{W}, j::_MM{W}) where {W} = vadd(i, vrange(j))
@inline Base.:(*)(i::_MM{W}, j::AbstractSIMDVector{W}) where {W} = vmul(vrange(i), j)
@inline Base.:(*)(i::AbstractSIMDVector{W}, j::_MM{W}) where {W} = vmul(i, vrange(j))


@inline vrange(::Val{W}) where {W} = vrange(Val{W}(), Float64)
@inline svrange(::Val{W}) where {W} = svrange(Val{W}(), Float64)

@inline vrange(i::_MM{W}, ::Type{Float64}) where {W} = vadd(vrange(Val{W}(), Float64), i.i)
@inline vrange(i::_MM{W}, ::Type{Float32}) where {W} = vadd(vrange(Val{W}(), Float32), i.i % Int32)
@inline vrange(i::_MM{W}, ::Type{Float16}) where {W} = vadd(vrange(Val{W}(), Float16), i.i % Int16)
@inline vrange(i::_MM{W}, ::Type{I}) where {W,I<:Integer} = vadd(vrange(Val{W}(), I), i.i % I)
@inline svrange(i::_MM, ::Type{T}) where {T} = SVec(vrange(i, T))


          
