
# Floating point arithmetic functions

for op ∈ (
    :(+), #:(-),
    :abs, :floor, :ceil, :round,
    :inv,# :log10, :log2, :exp2, :sin, :cos,
    :sqrt, :trunc
    # :exp, :log
    )
    rename = VECTOR_SYMBOLS[op]
    @eval begin
        @inline $rename(s1::FloatingTypes) = @fastmath $op(s1)

        @vectordef $rename function Base.$op(v1) where {W,T<:FloatingTypes}
            llvmwrap(Val{$(QuoteNode(op))}(), extract_data(v1))
        end
    end
end

@inline vsign(s1::FloatingTypes) = sign(s1)
@vectordef vsign function Base.sign(v1) where {W,T<:FloatingTypes}
    vifelse(
        visequal(extract_data(v1), vbroadcast(Vec{W,T}, zero(T))),
        vbroadcast(Vec{W,T}, zero(T)),
        vcopysign(vbroadcast(Vec{W,T}, one(T)), extract_data(v1))
    )
end
# @inline function vsign(v1::SVec{W,T}) where {W,T<:FloatingTypes}
#     SVec(vsign(extract_data(v1)))
# end
# @inline function Base.sign(v1::SVec{W,T}) where {W,T<:FloatingTypes}
#     SVec(vsign(extract_data(v1)))
# end

for op ∈ (
        :(+), :(-), :(*), :(/), :(%),# :(^),
        :copysign#, :max, :min
    )
    rename = VECTOR_SYMBOLS[op]
    # exact / explicit version
    erename = Symbol(:e, rename)
    @eval begin
        @inline $rename(s1::FloatingTypes, s2::FloatingTypes) = @fastmath $op(s1, s2)

        @vectordef $rename function Base.$op(v1, v2) where {W,T <: FloatingTypes}
            llvmwrap(Val{$(QuoteNode(op))}(), extract_data(v1), extract_data(v2))
        end

        @inline $erename(s1::FloatingTypes, s2::FloatingTypes) = $op(s1, s2)

        @evectordef $erename function Base.$op(v1, v2) where {W,T <: FloatingTypes}
            # llvmwrap_notfast(Val{$(QuoteNode(op))}(), extract_data(v1), extract_data(v2))
            llvmwrap_notfast(Val{$(QuoteNode(op))}(), v1, v2)
        end

        @inline function Base.$op(v1::AbstractStructVec{W,T1}, v2::AbstractStructVec{W,T2}) where {W,T1,T2}
            T = promote_type(T1, T2)
            $op(vconvert(SVec{W,T}, v1), vconvert(SVec{W,T}, v2))
        end
        @inline function $rename(v1::V1, v2::V2) where {W,V1<:AbstractSIMDVector{W},V2<:AbstractSIMDVector{W}}
            V = promote_vtype(V1, V2)
            $rename(vconvert(V, v1), vconvert(V, v2))
        end
    end
end
# Julia 1.0 was bad at resovling dispatch (it liked to pick the wrong method)
# On Julia 1.0, this was less specific:
# ::SVec{W,Float64}, ::SVec{W,Float64}
# Than:
# ::Union{Vec{W,T},AbstractStructVec{W,T}}, ::Union{Vec{W,T},AbstractStructVec{W,T}} where {T <: Union{Float64,Float32,...}
#
# Note that SVec{W,Float64} <: AbstractStructVec{W,Float64}.
# We have to make the type fully concrete (e.g., SVec{8,Float64}) before
# to get Julia 1.0 and 1.1 to dispatch on it.
# This is definitely not needed on Julia 1.4, so I place the bound there. I have not actually tested Julia 1.2 and Julia 1.3.
@static if VERSION ≤ v"1.4"
    for T ∈ [Float32, Float64, Int16, Int32, Int64, UInt16, UInt32, UInt64]
        for op ∈ (
            :(+), :(-), :(*), :(/), :(%),# :(^),
            :copysign#, :max, :min
        )
            rename = VECTOR_SYMBOLS[op]
            for W ∈ 1:16
                @eval @inline $rename(v1::SVec{$W,$T},v2::SVec{$W,$T}) = SVec($rename(extract_data(v1), extract_data(v2)))
            end
        end
        for W ∈ 1:16
            U = VectorizationBase.mask_type(W)
            @eval @inline vadd(m::Mask{$W,$U}, v::SVec{$W,$T}) = vifelse(m, vadd(v, one(v)), v)
            @eval @inline vadd(v::SVec{$W,$T}, m::Mask{$W,$U}) = vifelse(m, vadd(v, one(v)), v)
            @eval @inline vsub(m::Mask{$W,$U}, v::SVec{$W,$T}) = vsub(vifelse(m, one(v), zero(v)), v)
            @eval @inline vsub(v::SVec{$W,$T}, m::Mask{$W,$U}) = vifelse(m, vsub(v, one(v)), v)
            @eval @inline vmul(m::Mask{$W,$U}, v::SVec{$W,$T}) = vifelse(m, v, vzero(SVec{$W,$T}))
            @eval @inline vmul(v::SVec{$W,$T}, m::Mask{$W,$U}) = vifelse(m, v, vzero(SVec{$W,$T}))
        end
    end
end
@inline vfdiv(v1::_Vec{W,Int32}, v2::_Vec{W,Int32}) where {W} = vfdiv(vconvert(_Vec{W,Float32}, v1), vconvert(_Vec{W,Float32}, v2))
@inline vfdiv(v1::_Vec{W,Int64}, v2::_Vec{W,Int64}) where {W} = vfdiv(vconvert(_Vec{W,Float64}, v1), vconvert(_Vec{W,Float64}, v2))
@inline function Base.:(/)(v1::SVec{W,I},v2::SVec{W,I}) where {W,I<:Integer}
    T = sizeequivalentfloat(I)
    SVec(vfdiv(vconvert(Vec{W,T}, v1), vconvert(Vec{W,T}, v2)))
end
@inline Base.:(/)(v::SVec{W,I}, s::Integer) where {W,I} = v / vbroadcast(SVec{W,I}, s % I)
@inline Base.:(/)(s::Integer, v::SVec{W,I}) where {W,I} = vbroadcast(SVec{W,I}, s % I) / v
@inline vmax(x::Number, y::Number) = Base.FastMath.max_fast(x,y)
@vectordef vmax function Base.max(v1, v2) where {W,T<:FloatingTypes}
    vifelse(vgreater(extract_data(v1), extract_data(v2)), extract_data(v1), extract_data(v2))
end
@inline vmin(x::Number, y::Number) = Base.FastMath.min_fast(x,y)
@vectordef vmin function Base.min(v1, v2) where {W,T<:FloatingTypes}
    vifelse(vless(extract_data(v1), extract_data(v2)), extract_data(v1), extract_data(v2))
end


@inline vflipsign(s1::FloatingTypes, s2::FloatingTypes) = flipsign(s1, s2)

@vectordef vflipsign function Base.flipsign(v1, v2) where {W,T<:FloatingTypes}
    vifelse(vsignbit(extract_data(v2)), vsub(extract_data(v1)), extract_data(v1))
end

@inline vcopysign(v1::_Vec{W,T}, v2::_Vec{W,U}) where {W,T,U<:Unsigned} = vcopysign(v1, vreinterpret(_Vec{W,T}, v2))

for op ∈ (:fma, :muladd)
# let op = :fma
    rename = VECTOR_SYMBOLS[op]
    # if op == :muladd
        # rename = Symbol(:e,rename)
    # end
    @eval begin

        @vectordef $rename function Base.$op(v1, v2, v3) where {W,T<:FloatingTypes}
            llvmwrap(Val{$(QuoteNode(op))}(), extract_data(v1), extract_data(v2), extract_data(v3))
        end
        @vectordef $rename function Base.$op(s1::T, v2, v3) where {W,T<:FloatingTypes}
            llvmwrap(Val{$(QuoteNode(op))}(), vbroadcast(Vec{W,T}, s1), extract_data(v2), extract_data(v3))
        end
        @vectordef $rename function Base.$op(v1, s2::T, v3) where {W,T<:FloatingTypes}
            llvmwrap(Val{$(QuoteNode(op))}(), extract_data(v1), vbroadcast(Vec{W,T}, s2), extract_data(v3))
        end
        @vectordef $rename function Base.$op(v1, v2, s3::T) where {W,T<:FloatingTypes}
            llvmwrap(Val{$(QuoteNode(op))}(), extract_data(v1), extract_data(v2), vbroadcast(Vec{W,T}, s3))
        end
        @vectordef $rename function Base.$op(s1::T, s2::T, v3) where {W,T<:FloatingTypes}
            llvmwrap(Val{$(QuoteNode(op))}(), vbroadcast(Vec{W,T}, s1), vbroadcast(Vec{W,T}, s2), extract_data(v3))
        end
        @vectordef $rename function Base.$op(s1::T, v2, s3::T) where {W,T<:FloatingTypes}
            llvmwrap(Val{$(QuoteNode(op))}(), vbroadcast(Vec{W,T}, s1), extract_data(v2), vbroadcast(Vec{W,T}, s3))
        end
        @vectordef $rename function Base.$op(v1, s2::T, s3::T) where {W,T<:FloatingTypes}
            llvmwrap(Val{$(QuoteNode(op))}(), extract_data(v1), vbroadcast(Vec{W,T}, s2), vbroadcast(Vec{W,T}, s3))
        end
        @vpromote $rename 3
    end
end

# Type promotion



# Promote scalars of all ScalarTypes to vectors of FloatingTypes, leaving the
# vector type unchanged

for op ∈ (
        :(==), :(!=), :(<), :(<=), :(>), :(>=),
        :+, :-, :*, :/,# :^,
        :copysign, :flipsign, :max, :min, :%
    )
    rename = VECTOR_SYMBOLS[op]
    @eval begin
        @vectordef $rename function Base.$op(s1::ScalarTypes, v2) where {W,T<:FloatingTypes}
            ev2 = extract_data(v2)
            $rename(vbroadcast(typeof(ev2), s1), ev2)
        end
        @vectordef $rename function Base.$op(v1, s2::ScalarTypes) where {W,T<:FloatingTypes}
            ev1 = extract_data(v1)
            $rename(ev1, vbroadcast(typeof(ev1), s2))
        end
        @inline function $rename(s1::T2, v2::_Vec{W,T}) where {W,T<:Integer,T2<:FloatingTypes}
            $rename(vbroadcast(_Vec{W,T2}, s1), v2)
        end
        @inline function $rename(v1::_Vec{W,T}, s2::T2) where {W,T<:Integer,T2<:FloatingTypes}
            $rename(v1, vbroadcast(_Vec{W,T2}, s2))
        end
        @vpromote $rename 2
    end
end
for op ∈ (:(==), :(!=), :(<), :(<=), :(>), :(>=))
    rename = VECTOR_SYMBOLS[op]
    @eval begin
        @inline function Base.$op(s1::T2, v2::SVec{W,T}) where {W,T<:Integer,T2<:FloatingTypes}
            Mask{W}($rename(vbroadcast(Vec{W,T2}, s1), vconvert(Vec{W,T2}, extract_data(v2))))
        end
        @inline function Base.$op(v1::SVec{W,T}, s2::T2) where {W,T<:Integer,T2<:FloatingTypes}
            Mask{W}($rename(vconvert(Vec{W,T2}, extract_data(v1)), vbroadcast(Vec{W,T2}, s2)))
        end
        @inline function $rename(s1::T2, v2::SVec{W,T}) where {W,T<:Integer,T2<:FloatingTypes}
            Mask{W}($rename(vbroadcast(Vec{W,T2}, s1), extract_data(v2)))
        end
        @inline function $rename(v1::SVec{W,T}, s2::T2) where {W,T<:Integer,T2<:FloatingTypes}
            Mask{W}($rename(extract_data(v1), vbroadcast(Vec{W,T2}, s2)))
        end
    end
end
for op ∈ (:+, :-, :*, :/, :copysign, :flipsign, :max, :min, :%)
    rename = VECTOR_SYMBOLS[op]
    @eval begin
        @inline function Base.$op(s1::T2, v2::SVec{W,T}) where {W,T<:Integer,T2<:FloatingTypes}
            SVec($rename(vbroadcast(Vec{W,T2}, s1), vconvert(Vec{W,T2}, extract_data(v2))))
        end
        @inline function Base.$op(v1::SVec{W,T}, s2::T2) where {W,T<:Integer,T2<:FloatingTypes}
            SVec($rename(vconvert(Vec{W,T2}, extract_data(v1)), vbroadcast(Vec{W,T2}, s2)))
        end
        @inline function $rename(s1::T2, v2::SVec{W,T}) where {W,T<:Integer,T2<:FloatingTypes}
            SVec($rename(vbroadcast(Vec{W,T2}, s1), extract_data(v2)))
        end
        @inline function $rename(v1::SVec{W,T}, s2::T2) where {W,T<:Integer,T2<:FloatingTypes}
            SVec($rename(extract_data(v1), vbroadcast(Vec{W,T2}, s2)))
        end
    end
end

for op ∈ (
        :(+), :(-), :(*), :(/), :(%),# :(^),
        :copysign#, :max, :min
    )
    # exact / explicit version
    erename = Symbol(:e, VECTOR_SYMBOLS[op])
    @eval begin

        @evectordef $erename function Base.$op(v1, s2::ScalarTypes) where {W,T <: FloatingTypes}
            $erename(extract_data(v1), vbroadcast(Vec{W,T}, s2))
        end
        @evectordef $erename function Base.$op(s1::ScalarTypes, v2) where {W,T <: FloatingTypes}
            $erename(vbroadcast(Vec{W,T}, s1), extract_data(v2))
        end
    end
end
 
# Reduction operations

# TODO: map, mapreduce

function getneutral(op::Symbol, ::Type{T}) where {T}
    if op == :&
        return ~zero(T)
    elseif op == :|
        return zero(T)
    elseif op == :max
        return typemin(T)
    elseif op == :min
        return typemax(T)
    elseif op == :+
        return zero(T)
    elseif op == :*
        return one(T)
    end
    throw("Op $op not recognized.")
end


@static if Base.libllvm_version >= v"9"
    @generated function vsum(v::_Vec{_W,T}) where {_W,T<:FloatingTypes}
        W = _W + 1
        instrs = String[]
        typ = llvmtype(T)
        vtyp = "<$W x $typ>"
        bits = 8sizeof(T)
        ins = "@llvm.experimental.vector.reduce.v2.fadd.f$(bits).v$(W)f$(bits)"
        decl = "declare $(typ) $(ins)($(typ), $(vtyp))"
        push!(instrs, "%res = call $(fastflags(T)) $typ $ins($typ 0.0, $vtyp %0)")
        push!(instrs, "ret $typ %res")
        quote
            $(Expr(:meta, :inline))
            llvmcall(
                $((decl, join(instrs, "\n"))),
                $T, Tuple{Vec{$W,$T}}, v
            )
        end
    end
    @generated function vprod(v::_Vec{_W,T}) where {_W,T<:FloatingTypes}
        W = _W + 1
        instrs = String[]
        typ = llvmtype(T)
        vtyp = "<$W x $typ>"
        bits = 8sizeof(T)
        ins = "@llvm.experimental.vector.reduce.v2.fmul.f$(bits).v$(W)f$(bits)"
        decl = "declare $(typ) $(ins)($(typ), $(vtyp))"
        push!(instrs, "%res = call $(fastflags(T)) $typ $ins($typ 1.0, $vtyp %0)")
        push!(instrs, "ret $typ %res")
        quote
            $(Expr(:meta, :inline))
            llvmcall(
                $((decl, join(instrs, "\n"))),
                $T, Tuple{Vec{$W,$T}}, v
            )
        end
    end
    @generated function reduced_add(v::_Vec{_W,T}, s::T) where {_W,T<:FloatingTypes}
        W = _W + 1
        instrs = String[]
        typ = llvmtype(T)
        vtyp = "<$W x $typ>"
        bits = 8sizeof(T)
        ins = "@llvm.experimental.vector.reduce.v2.fadd.f$(bits).v$(W)f$(bits)"
        decl = "declare $(typ) $(ins)($(typ), $(vtyp))"
        push!(instrs, "%res = call $(fastflags(T)) $typ $ins($typ %1, $vtyp %0)")
        push!(instrs, "ret $typ %res")
        quote
            $(Expr(:meta, :inline))
            llvmcall(
                $((decl, join(instrs, "\n"))),
                $T, Tuple{Vec{$W,$T},$T}, v, s
            )
        end
    end
    @generated function reduced_prod(v::_Vec{_W,T}, s::T) where {_W,T<:FloatingTypes}
        W = _W + 1
        instrs = String[]
        typ = llvmtype(T)
        vtyp = "<$W x $typ>"
        bits = 8sizeof(T)
        ins = "@llvm.experimental.vector.reduce.v2.fmul.f$(bits).v$(W)f$(bits)"
        decl = "declare $(typ) $(ins)($(typ), $(vtyp))"
        push!(instrs, "%res = call $(fastflags(T)) $typ $ins($typ %1, $vtyp %0)")
        push!(instrs, "ret $typ %res")
        quote
            $(Expr(:meta, :inline))
            llvmcall(
                $((decl, join(instrs, "\n"))),
                $T, Tuple{Vec{$W,$T},$T}, v, s
            )
        end
    end
    @generated function vsum(v::_Vec{_W,T}) where {_W,T<:IntegerTypes}
        W = _W + 1
        instrs = String[]
        typ = llvmtype(T)
        vtyp = "<$W x $typ>"
        bits = 8sizeof(T)
        ins = "@llvm.experimental.vector.reduce.add.v$(W)i$(bits)"
        decl = "declare $(typ) $(ins)($(vtyp))"
        push!(instrs, "%res = call $typ $ins($vtyp %0)")
        push!(instrs, "ret $typ %res")
        quote
            $(Expr(:meta, :inline))
            llvmcall(
                $((decl, join(instrs, "\n"))),
                $T, Tuple{Vec{$W,$T}}, v
            )
        end
    end
    @generated function vprod(v::_Vec{_W,T}) where {_W,T<:IntegerTypes}
        W = _W + 1
        instrs = String[]
        typ = llvmtype(T)
        vtyp = "<$W x $typ>"
        bits = 8sizeof(T)
        ins = "@llvm.experimental.vector.reduce.mul.v$(W)i$(bits)"
        decl = "declare $(typ) $(ins)($(vtyp))"
        push!(instrs, "%res = call $typ $ins($vtyp %0)")
        push!(instrs, "ret $typ %res")
        quote
            $(Expr(:meta, :inline))
            llvmcall(
                $((decl, join(instrs, "\n"))),
                $T, Tuple{Vec{$W,$T}}, v
            )
        end
    end
else
    @generated function llvmwrapreduce(::Val{Op}, v::_Vec{_W,T}) where {Op,_W,T}
        @assert isa(Op, Symbol)
        W = _W + 1
        z = getneutral(Op, T)
        typ = llvmtype(T)
        decls = String[]
        instrs = String[]
        n = W
        nam = "%0"
        nold,n = n, VectorizationBase.nextpow2(n)
        if n > nold
            namold,nam = nam,"%vec_$n"
            append!(instrs,
                    extendvector(namold, nold, typ, n, n-nold,
                                 llvmtypedconst(T,z), nam))
        end
        while n > 1
            nold,n = n, div(n, 2)
            namold,nam = nam,"%vec_$n"
            vtyp = "<$n x $typ>"
            ins = llvmins(Op, n, T)
            append!(instrs, subvector(namold, nold, typ, "$(nam)_1", n, 0))
            append!(instrs, subvector(namold, nold, typ, "$(nam)_2", n, n))
            if ins[1] == '@'
                push!(decls, "declare $vtyp $ins($vtyp, $vtyp)")
                push!(instrs,
                      "$nam = call $vtyp $ins($vtyp $(nam)_1, $vtyp $(nam)_2)")
            else
                push!(instrs, "$nam = $ins $vtyp $(nam)_1, $(nam)_2")
            end
        end
        push!(instrs, "%res = extractelement <$n x $typ> $nam, i32 0")
        push!(instrs, "ret $typ %res")
        quote
            $(Expr(:meta, :inline))
            llvmcall(
                $((join(decls, "\n"), join(instrs, "\n"))),
                $T, Tuple{Vec{$W,$T}}, v
            )
        end
    end

    
    for (name, rename, op) ∈ [(:(Base.sum),:vsum,:+), (:(Base.prod),:vprod,:*)]
        @eval begin
            @inline $rename(v::AbstractSIMDVector{W}) where {W} = llvmwrapreduce(Val{$(QuoteNode(op))}(), extract_data(v))
        end
    end

end
@static if Base.libllvm_version >= v"8"
    
    @generated function vsub(v::_Vec{_W,T}) where {_W,T<:FloatingTypes}
        W = _W + 1
        typ = llvmtype(T)
        vtyp = "<$W x $typ>"
        instrs = "%res = fneg $(fastflags(T)) $vtyp %0\nret $vtyp %res"
        quote
            $(Expr(:meta, :inline))
            llvmcall( $instrs, Vec{$W,$T}, Tuple{Vec{$W,$T}}, v )
        end
    end
    @inline Base.:(-)(v::SVec{W,T}) where {W,T} = SVec(vsub(extract_data(v)))
    @inline vsub(v::SVec{W,T}) where {W,T} = SVec(vsub(extract_data(v)))
else
    let rename = VECTOR_SYMBOLS[:-]
        @eval begin
            @vectordef $rename function Base.:(-)(v1) where {W,T<:FloatingTypes}
                llvmwrap(Val{:(-)}(), extract_data(v1))
            end
        end
    end

end
@inline vsub(x::FloatingTypes) = Base.FastMath.sub_fast(x)

@static if Base.libllvm_version >= v"6"
    @generated function vmaximum(v::_Vec{_W,T}) where {_W,T<:ScalarTypes}
        W = _W + 1
        instrs = String[]
        typ = llvmtype(T)
        vtyp = "<$W x $typ>"
        bits = 8sizeof(T)
        if T <: FloatingTypes
            prefix = prefix2 = 'f'
        else
            prefix2 = 'i'
            prefix = T <: Signed ? 's' : 'u'
        end
        ins = "@llvm.experimental.vector.reduce.$(prefix)max.v$(W)$(prefix2)$(bits)"
        decl = "declare $(typ) $(ins)($(vtyp))"
        push!(instrs, "%res = call $typ $ins($vtyp %0)")
        push!(instrs, "ret $typ %res")
        quote
            $(Expr(:meta, :inline))
            llvmcall(
                $((decl, join(instrs, "\n"))),
                $T, Tuple{Vec{$W,$T}}, v
            )
        end
    end
    @generated function vminimum(v::_Vec{_W,T}) where {_W,T<:ScalarTypes}
        W = _W + 1
        instrs = String[]
        typ = llvmtype(T)
        vtyp = "<$W x $typ>"
        bits = 8sizeof(T)
        if T <: FloatingTypes
            prefix = prefix2 = 'f'
        else
            prefix2 = 'i'
            prefix = T <: Signed ? 's' : 'u'
        end
        ins = "@llvm.experimental.vector.reduce.$(prefix)min.v$(W)$(prefix2)$(bits)"
        decl = "declare $(typ) $(ins)($(vtyp))"
        push!(instrs, "%res = call $typ $ins($vtyp %0)")
        push!(instrs, "ret $typ %res")
        quote
            $(Expr(:meta, :inline))
            llvmcall(
                $((decl, join(instrs, "\n"))),
                $T, Tuple{Vec{$W,$T}}, v
            )
        end
    end
else
    for (name, rename, op) ∈ [(:(Base.maximum), :vmaximum, :max), (:(Base.minimum), :vminimum, :min)]
                              
        @eval begin
            @inline $rename(v::AbstractSIMDVector{W}) where {W} = llvmwrapreduce(Val{$(QuoteNode(op))}(), extract_data(v))
        end
    end

end
for (name, rename, op) ∈ [(:(Base.maximum), :vmaximum, :max), (:(Base.minimum), :vminimum, :min), (:(Base.sum),:vsum,:+), (:(Base.prod),:vprod,:*)]
    @eval begin
        @inline $name(v::SVec{W,T}) where {W,T} = $rename(extract_data(v))
        @inline $rename(v::SVec{W,T}) where {W,T} = $rename(extract_data(v))
        # @inline $name(v::SVec{W,T}) where {W,T} = llvmwrapreduce(Val{$(QuoteNode(op))}(), extract_data(v))
    end
end



# TODO: Handle cases with vectors of different lengths correctly!
@inline vmul(x::Number, y::Number) = Base.FastMath.mul_fast(x, y)
@inline vadd(x::Number, y::Number) = Base.FastMath.add_fast(x, y)
@inline vmul(x,y,z...) = vmul(x,vmul(y,z...))
@inline vadd(x,y,z...) = vadd(x,vadd(y,z...))
@inline vadd(a,b,c,d) = evadd(vadd(a,b), vadd(c,d))
@inline vadd(a,b,c,d,e,f) = evadd(evadd(vadd(a,b), vadd(c,d)),vad(e,f))
@inline vadd(a,b,c,d,e,f,g,h) = vadd(evadd(vadd(a,b), vadd(c,d)),evadd(vadd(e,f),vadd(g,h)))
@inline Base.:(+)(a::SVec, b, c, d) = evadd(vadd(a,b), vadd(c,d))
@inline vmuladd(a::Number, b::Number, c::Number) = muladd(a, b, c)
# These intrinsics are not FastMath.

@generated function vfmadd(v1::_Vec{_W,T}, v2::_Vec{_W,T}, v3::_Vec{_W,T}) where {_W,T<:FloatingTypes}
    W = _W + 1
    typ = llvmtype(T)
    vtyp = "<$W x $typ>"
    # ins = "@llvm.fma.v$(W)f$(8*sizeof(T))"
    ins = "@llvm.fmuladd.v$(W)f$(8*sizeof(T))"
    decl = "declare $vtyp $ins($vtyp, $vtyp, $vtyp)"
    instrs = "%res = call $vtyp $ins($vtyp %0, $vtyp %1, $vtyp %2)\nret $vtyp %res"
    quote
        $(Expr(:meta, :inline))
        llvmcall(
            $((decl,instrs)),
            Vec{$W,$T}, Tuple{Vec{$W,$T}, Vec{$W,$T}, Vec{$W,$T}},
            v1, v2, v3
        )
    end
end
vfmadd(a::Number, b::Number, c::Number) = muladd(a, b, c)
vfnmadd(a::Number, b::Number, c::Number) = muladd(-a, b, c)
vfmsub(a::Number, b::Number, c::Number) = muladd(a, b, -c)
vfnmsub(a::Number, b::Number, c::Number) = -muladd(a, b, c)
@inline vfmadd(a::_Vec{W,T}, b::_Vec{W,T}, c::_Vec{W,T}) where {W,T <: Integer} = vmuladd(a, b, c)
@inline vfnmadd(a::_Vec{W,T}, b::_Vec{W,T}, c::_Vec{W,T}) where {W,T} = vfmadd(vsub(a), b, c)
@inline vfmsub(a::_Vec{W,T}, b::_Vec{W,T}, c::_Vec{W,T}) where {W,T} = vfmadd(a, b, vsub(c))
@inline vfnmsub(a::_Vec{W,T}, b::_Vec{W,T}, c::_Vec{W,T}) where {W,T} = vsub(vfmadd(a, b, c))
@vpromote vfmadd 3
@vpromote vfnmadd 3
@vpromote vfmsub 3
@vpromote vfnmsub 3

@generated function vfmadd_fast(v1::_Vec{_W,T}, v2::_Vec{_W,T}, v3::_Vec{_W,T}) where {_W,T<:FloatingTypes}
    W = _W + 1
    typ = llvmtype(T)
    vtyp = "<$W x $typ>"
    # ins = "@llvm.fmuladd.v$(W)f$(8*sizeof(T))"
    # decls = "declare $vtyp $ins($vtyp, $vtyp, $vtyp)"
    # instrs = "%res = call fast $vtyp $ins($vtyp %0, $vtyp %1, $vtyp %2)\nret $vtyp %res"
    # I don't really want reassoc on the fmul part, so this is what I'm going with.
    instrs = """
    %prod = fmul nnan ninf nsz arcp contract $vtyp %0, %1
    %res = fadd $(fastflags(T)) $vtyp %prod, %2
    ret $vtyp %res
    """
    quote
        $(Expr(:meta, :inline))
        llvmcall(
            $instrs, Vec{$W,$T}, Tuple{Vec{$W,$T}, Vec{$W,$T}, Vec{$W,$T}}, v1, v2, v3
        )
    end
end
vfmadd_fast(a::Number, b::Number, c::Number) = Base.FastMath.add_fast(Base.FastMath.mul_fast(a, b), c)
vfnmadd_fast(a::Number, b::Number, c::Number) = Base.FastMath.sub_fast(c, Base.FastMath.mul_fast(a, b))
vfmsub_fast(a::Number, b::Number, c::Number) = Base.FastMath.sub_fast(Base.FastMath.mul_fast(a, b), c)
vfnmsub_fast(a::Number, b::Number, c::Number) = Base.FastMath.sub_fast(Base.FastMath.add_fast(Base.FastMath.mul_fast(a, b), c))
@inline vfmadd_fast(a::_Vec{W,T}, b::_Vec{W,T}, c::_Vec{W,T}) where {W,T <: Integer} = vmuladd(a, b, c)
@inline vfnmadd_fast(a::_Vec{W,T}, b::_Vec{W,T}, c::_Vec{W,T}) where {W,T} = vfmadd_fast(vsub(a), b, c)
@inline vfmsub_fast(a::_Vec{W,T}, b::_Vec{W,T}, c::_Vec{W,T}) where {W,T} = vfmadd_fast(a, b, vsub(c))
@inline vfnmsub_fast(a::_Vec{W,T}, b::_Vec{W,T}, c::_Vec{W,T}) where {W,T} = vsub(vfmadd_fast(a, b, c))
@inline vfmadd_fast(m::AbstractMask{W}, b::AbstractSIMDVector{W}, c::AbstractSIMDVector{W}) where {W} = vifelse(m, vadd(b, c), c)
@inline vfmadd_fast(b::AbstractSIMDVector{W}, m::AbstractMask{W}, c::AbstractSIMDVector{W}) where {W} = vifelse(m, vadd(b, c), c)
@inline vfnmadd_fast(m::AbstractMask{W}, b::AbstractSIMDVector{W}, c::AbstractSIMDVector{W}) where {W} = vifelse(m, vsub(c, b), c)
@inline vfnmadd_fast(b::AbstractSIMDVector{W}, m::AbstractMask{W}, c::AbstractSIMDVector{W}) where {W} = vifelse(m, vsub(c, b), c)
@inline vfmsub_fast(m::AbstractMask{W}, b::V, c::AbstractSIMDVector{W}) where {W,V<:AbstractSIMDVector{W}} = vsub(vifelse(m, b, vzero(V)), c)
@inline vfmsub_fast(b::V, m::AbstractMask{W}, c::AbstractSIMDVector{W}) where {W,V<:AbstractSIMDVector{W}} = vsub(vifelse(m, b, vzero(V)), c)
@inline vfnmsub_fast(m::AbstractMask{W}, b::AbstractSIMDVector{W}, c::AbstractSIMDVector{W}) where {W} = vsub(vfmadd_fast(m, b, c))
@inline vfnmsub_fast(b::AbstractSIMDVector{W}, m::AbstractMask{W}, c::AbstractSIMDVector{W}) where {W} = vsub(vfmadd_fast(m, b, c))
@inline vfmadd(m::AbstractMask{W}, b::AbstractSIMDVector{W}, c::AbstractSIMDVector{W}) where {W} = vifelse(m, vadd(b, c), c)
@inline vfmadd(b::AbstractSIMDVector{W}, m::AbstractMask{W}, c::AbstractSIMDVector{W}) where {W} = vifelse(m, vadd(b, c), c)
@inline vfnmadd(m::AbstractMask{W}, b::AbstractSIMDVector{W}, c::AbstractSIMDVector{W}) where {W} = vifelse(m, vsub(c, b), c)
@inline vfnmadd(b::AbstractSIMDVector{W}, m::AbstractMask{W}, c::AbstractSIMDVector{W}) where {W} = vifelse(m, vsub(c, b), c)
@inline vfmsub(m::AbstractMask{W}, b::V, c::AbstractSIMDVector{W}) where {W,V<:AbstractSIMDVector{W}} = vsub(vifelse(m, b, vzero(V)), c)
@inline vfmsub(b::V, m::AbstractMask{W}, c::AbstractSIMDVector{W}) where {W,V<:AbstractSIMDVector{W}} = vsub(vifelse(m, b, vzero(V)), c)
@inline vfnmsub(m::AbstractMask{W}, b::AbstractSIMDVector{W}, c::AbstractSIMDVector{W}) where {W} = vsub(vfmadd(m, b, c))
@inline vfnmsub(b::AbstractSIMDVector{W}, m::AbstractMask{W}, c::AbstractSIMDVector{W}) where {W} = vsub(vfmadd(m, b, c))
@inline vfmadd231(m::AbstractMask{W}, b::AbstractSIMDVector{W}, c::AbstractSIMDVector{W}) where {W} = vifelse(m, vadd(b, c), c)
@inline vfmadd231(b::AbstractSIMDVector{W}, m::AbstractMask{W}, c::AbstractSIMDVector{W}) where {W} = vifelse(m, vadd(b, c), c)
@inline vfnmadd231(m::AbstractMask{W}, b::AbstractSIMDVector{W}, c::AbstractSIMDVector{W}) where {W} = vifelse(m, vsub(c, b), c)
@inline vfnmadd231(b::AbstractSIMDVector{W}, m::AbstractMask{W}, c::AbstractSIMDVector{W}) where {W} = vifelse(m, vsub(c, b), c)
@inline vfmsub231(m::AbstractMask{W}, b::V, c::AbstractSIMDVector{W}) where {W,V<:AbstractSIMDVector{W}} = vsub(vifelse(m, b, vzero(V)), c)
@inline vfmsub231(b::V, m::AbstractMask{W}, c::AbstractSIMDVector{W}) where {W,V<:AbstractSIMDVector{W}} = vsub(vifelse(m, b, vzero(V)), c)
@inline vfnmsub231(m::AbstractMask{W}, b::AbstractSIMDVector{W}, c::AbstractSIMDVector{W}) where {W} = vsub(vfmadd231(m, b, c))
@inline vfnmsub231(b::AbstractSIMDVector{W}, m::AbstractMask{W}, c::AbstractSIMDVector{W}) where {W} = vsub(vfmadd231(m, b, c))
@vpromote vfmadd_fast 3
@vpromote vfnmadd_fast 3
@vpromote vfmsub_fast 3
@vpromote vfnmsub_fast 3

@inline Base.:(+)(m::AbstractMask{W}, v::AbstractStructVec{W}) where {W} = vifelse(m, vadd(v, one(v)), v)
@inline Base.:(+)(v::AbstractStructVec{W}, m::AbstractMask{W}) where {W} = vifelse(m, vadd(v, one(v)), v)
@inline vadd(m::AbstractMask{W}, v::AbstractSIMDVector{W}) where {W} = vifelse(m, vadd(v, one(v)), v)
@inline vadd(v::AbstractSIMDVector{W}, m::AbstractMask{W}) where {W} = vifelse(m, vadd(v, one(v)), v)
@inline Base.:(-)(m::AbstractMask{W}, v::AbstractStructVec{W}) where {W} = vsub(vifelse(m, one(v), zero(v)), v)
@inline Base.:(-)(v::AbstractStructVec{W}, m::AbstractMask{W}) where {W} = vifelse(m, vsub(v, one(v)), v)
@inline vsub(m::AbstractMask{W}, v::AbstractSIMDVector{W}) where {W} = vsub(vifelse(m, one(v), zero(v)), v)
@inline vsub(v::AbstractSIMDVector{W}, m::AbstractMask{W}) where {W} = vifelse(m, vsub(v, one(v)), v)
@inline Base.:(*)(v::AbstractStructVec{W,T}, m::Mask{W}) where {W,T} = SVec(vifelse(m, extract_data(v), vzero(Vec{W,T})))
@inline Base.:(*)(m::Mask{W}, v::AbstractStructVec{W,T}) where {W,T} = SVec(vifelse(m, extract_data(v), vzero(Vec{W,T})))
@inline vmul(v::AbstractStructVec{W,T}, m::Mask{W}) where {W,T} = SVec(vifelse(m, extract_data(v), vzero(Vec{W,T})))
@inline vmul(m::Mask{W}, v::AbstractStructVec{W,T}) where {W,T} = SVec(vifelse(m, extract_data(v), vzero(Vec{W,T})))
@inline Base.:(*)(v::AbstractStructVec{W,T}, m::SVec{W,Bool}) where {W,T} = SVec(vifelse(extract_data(m), extract_data(v), vzero(Vec{W,T})))
@inline Base.:(*)(m::SVec{W,Bool}, v::AbstractStructVec{W,T}) where {W,T} = SVec(vifelse(extract_data(m), extract_data(v), vzero(Vec{W,T})))
@inline vmul(v::AbstractStructVec{W,T}, m::SVec{W,Bool}) where {W,T} = SVec(vifelse(extract_data(m), extract_data(v), vzero(Vec{W,T})))
@inline vmul(m::SVec{W,Bool}, v::AbstractStructVec{W,T}) where {W,T} = SVec(vifelse(extract_data(m), extract_data(v), vzero(Vec{W,T})))

# Lowers to same split mul-add llvm as the 
# @inline vfmadd(a, b, c) = vadd(vmul( a, b), c)
# definition, so I wont bother implementing
# vfnmadd, vfmsub, and vfnmsub
# in this manner.

# @inline Base.:*(a::IntegerTypes, b::SVec{W,T}) where {W,T} = SVec{W,T}(a) * b
# @inline Base.:*(a::T, b::SVec{W,<:IntegerTypes}) where {W,T<:FloatingTypes} = SVec(vmul(vbroadcast(Vec{W,T}, a), vconvert(Vec{W,T}, extract_data(b))))

# const Vec{W,T} = NTuple{W,Core.VecElement{T}}
@inline vfmadd231(a::_Vec, b, c) = vfmadd(a, b, c)
@inline vfnmadd231(a::_Vec, b, c) = vfnmadd(a, b, c)
@inline vfmsub231(a::_Vec, b, c) = vfmsub(a, b, c)
@inline vfnmsub231(a::_Vec, b, c) = vfnmsub(a, b, c)
@inline vfmadd231(a::Number, b::Number, c::Number) = Base.FastMath.add_fast(Base.FastMath.mul_fast(a, b), c)
@inline vfnmadd231(a::Number, b::Number, c::Number) = Base.FastMath.sub_fast(c, Base.FastMath.mul_fast(a, b))
@inline vfmsub231(a::Number, b::Number, c::Number) = Base.FastMath.sub_fast(Base.FastMath.mul_fast(a, b), c)
@inline vfnmsub231(a::Number, b::Number, c::Number) = Base.FastMath.sub_fast(Base.FastMath.sub_fast(c), Base.FastMath.mul_fast(a, b))
@inline function vifelse(::typeof(vfmadd231), m::AbstractMask, a, b, c)
    vifelse(m, vfmadd(a, b, c), c)
end
@inline function vifelse(::typeof(vfnmadd231), m::AbstractMask, a, b, c)
    vifelse(m, vfnmadd(a, b, c), c)
end
@inline function vifelse(::typeof(vfmsub231), m::AbstractMask, a, b, c)
    vifelse(m, vfmsub(a, b, c), c)
end
@inline function vifelse(::typeof(vfnmsub231), m::AbstractMask, a, b, c)
    vifelse(m, vfnmsub(a, b, c), c)
end
@inline function vifelse(::typeof(vfmadd_fast), m::AbstractMask, a, b, c)
    vifelse(m, vfmadd(a, b, c), c)
end
@inline function vifelse(::typeof(vfnmadd_fast), m::AbstractMask, a, b, c)
    vifelse(m, vfnmadd(a, b, c), c)
end
@inline function vifelse(::typeof(vfmsub_fast), m::AbstractMask, a, b, c)
    vifelse(m, vfmsub(a, b, c), c)
end
@inline function vifelse(::typeof(vfnmsub_fast), m::AbstractMask, a, b, c)
    vifelse(m, vfnmsub(a, b, c), c)
end

@vpromote vfmadd231 3
@vpromote vfnmadd231 3
@vpromote vfmsub231 3
@vpromote vfnmsub231 3

if VectorizationBase.FMA3
    for T ∈ [Float32,Float64]
        W = 16 ÷ sizeof(T)
        local suffix = T == Float32 ? "ps" : "pd"
        typ = llvmtype(T)
        while W <= VectorizationBase.REGISTER_SIZE ÷ sizeof(T)
            vfmadd_str = """%res = call <$W x $(typ)> asm "vfmadd231$(suffix) \$3, \$2, \$1", "=v,0,v,v"(<$W x $(typ)> %2, <$W x $(typ)> %1, <$W x $(typ)> %0)
                ret <$W x $(typ)> %res"""
            vfnmadd_str = """%res = call <$W x $(typ)> asm "vfnmadd231$(suffix) \$3, \$2, \$1", "=v,0,v,v"(<$W x $(typ)> %2, <$W x $(typ)> %1, <$W x $(typ)> %0)
                ret <$W x $(typ)> %res"""
            vfmsub_str = """%res = call <$W x $(typ)> asm "vfmsub231$(suffix) \$3, \$2, \$1", "=v,0,v,v"(<$W x $(typ)> %2, <$W x $(typ)> %1, <$W x $(typ)> %0)
                ret <$W x $(typ)> %res"""
            vfnmsub_str = """%res = call <$W x $(typ)> asm "vfnmsub231$(suffix) \$3, \$2, \$1", "=v,0,v,v"(<$W x $(typ)> %2, <$W x $(typ)> %1, <$W x $(typ)> %0)
                ret <$W x $(typ)> %res"""
            @eval begin
                @inline function vfmadd231(a::Vec{$W,$T}, b::Vec{$W,$T}, c::Vec{$W,$T})
                    llvmcall($vfmadd_str, Vec{$W,$T}, Tuple{Vec{$W,$T},Vec{$W,$T},Vec{$W,$T}}, a, b, c)
                end
                @inline function vfnmadd231(a::Vec{$W,$T}, b::Vec{$W,$T}, c::Vec{$W,$T})
                    llvmcall($vfnmadd_str, Vec{$W,$T}, Tuple{Vec{$W,$T},Vec{$W,$T},Vec{$W,$T}}, a, b, c)
                end
                @inline function vfmsub231(a::Vec{$W,$T}, b::Vec{$W,$T}, c::Vec{$W,$T})
                    llvmcall($vfmsub_str, Vec{$W,$T}, Tuple{Vec{$W,$T},Vec{$W,$T},Vec{$W,$T}}, a, b, c)
                end
                @inline function vfnmsub231(a::Vec{$W,$T}, b::Vec{$W,$T}, c::Vec{$W,$T})
                    llvmcall($vfnmsub_str, Vec{$W,$T}, Tuple{Vec{$W,$T},Vec{$W,$T},Vec{$W,$T}}, a, b, c)
                end
            end
            if VectorizationBase.AVX512BW && W ≥ 8
                vfmaddmask_str = """%res = call <$W x $(typ)> asm "vfmadd231$(suffix) \$3, \$2, \$1 {\$4}", "=v,0,v,v,^Yk"(<$W x $(typ)> %2, <$W x $(typ)> %1, <$W x $(typ)> %0, i$W %3)
                    ret <$W x $(typ)> %res"""
                vfnmaddmask_str = """%res = call <$W x $(typ)> asm "vfnmadd231$(suffix) \$3, \$2, \$1 {\$4}", "=v,0,v,v,^Yk"(<$W x $(typ)> %2, <$W x $(typ)> %1, <$W x $(typ)> %0, i$W %3)
                    ret <$W x $(typ)> %res"""
                vfmsubmask_str = """%res = call <$W x $(typ)> asm "vfmsub231$(suffix) \$3, \$2, \$1 {\$4}", "=v,0,v,v,^Yk"(<$W x $(typ)> %2, <$W x $(typ)> %1, <$W x $(typ)> %0, i$W %3)
                    ret <$W x $(typ)> %res"""
                vfnmsubmask_str = """%res = call <$W x $(typ)> asm "vfnmsub231$(suffix) \$3, \$2, \$1 {\$4}", "=v,0,v,v,^Yk"(<$W x $(typ)> %2, <$W x $(typ)> %1, <$W x $(typ)> %0, i$W %3)
                    ret <$W x $(typ)> %res"""
                U = VectorizationBase.mask_type(W)
                @eval begin
                    @inline function vifelse(::typeof(vfmadd231), m::Mask{$W,$U}, a::SVec{$W,$T}, b::SVec{$W,$T}, c::SVec{$W,$T})
                        SVec(llvmcall($vfmaddmask_str, Vec{$W,$T}, Tuple{Vec{$W,$T},Vec{$W,$T},Vec{$W,$T},$U}, extract_data(a), extract_data(b), extract_data(c), extract_data(m)))
                    end
                    @inline function vifelse(::typeof(vfnmadd231), m::Mask{$W,$U}, a::SVec{$W,$T}, b::SVec{$W,$T}, c::SVec{$W,$T})
                        SVec(llvmcall($vfnmaddmask_str, Vec{$W,$T}, Tuple{Vec{$W,$T},Vec{$W,$T},Vec{$W,$T},$U}, extract_data(a), extract_data(b), extract_data(c), extract_data(m)))
                    end
                    @inline function vifelse(::typeof(vfmsub231), m::Mask{$W,$U}, a::SVec{$W,$T}, b::SVec{$W,$T}, c::SVec{$W,$T})
                        SVec(llvmcall($vfmsubmask_str, Vec{$W,$T}, Tuple{Vec{$W,$T},Vec{$W,$T},Vec{$W,$T},$U}, extract_data(a), extract_data(b), extract_data(c), extract_data(m)))
                    end
                    @inline function vifelse(::typeof(vfnmsub231), m::Mask{$W,$U}, a::SVec{$W,$T}, b::SVec{$W,$T}, c::SVec{$W,$T})
                        SVec(llvmcall($vfnmsubmask_str, Vec{$W,$T}, Tuple{Vec{$W,$T},Vec{$W,$T},Vec{$W,$T},$U}, extract_data(a), extract_data(b), extract_data(c), extract_data(m)))
                    end
                end
            end
            W += W
        end
    end
end

# @generated function rsqrt_fast(x::NTuple{16,Core.VecElement{Float32}})
#     if VectorizationBase.REGISTER_SIZE == 64
#         return quote
#             $(Expr(:meta,:inline))
#             llvmcall(
#             """ %rs = call <16 x float> asm "vrsqrt14ps \$1, \$0", "=x,x"(<16 x float> %0)
#                 ret <16 x float> %rs""",
#             NTuple{16,Core.VecElement{Float32}}, Tuple{NTuple{16,Core.VecElement{Float32}}}, x)
#         end
#     else
#         return quote
#             vinv(vsqrt(x))
#         end
#     end
# end
# @inline function rsqrt(x::NTuple{16,Core.VecElement{Float32}})
#     r = rsqrt_fast(x)
#     # Performs a Newton step to increase accuracy.
#     # ns = vmuladd(vmul(-0.5f0, x), vmul(r, r), 1.5f0)
#     # vmul(r, ns)
#     ns = vfma(vmul(r,r), x, -3.0f0)
#     vmul(vmul(-0.5f0, r), ns)
# end
@inline rsqrt_fast(v) = vinv(vsqrt(v))
@inline rsqrt_fast(x::SVec) = SVec(rsqrt_fast(extract_data(x)))
@inline rsqrt(x::SVec) = SVec(rsqrt(extract_data(x)))
@inline rsqrt(x) = vinv(vsqrt(x))
@inline vinv(x::AbstractStructVec) = SVec(vinv(extract_data(x)))
@inline vinv(x::IntegerTypes) = vinv(float(x))
@inline vinv(x::_Vec{W,I}) where {W, I <: Union{Int64,UInt64}} = evfdiv(vone(_Vec{W,Float64}), vconvert(_Vec{W,Float64}, x))
@inline vinv(x::AbstractStructVec{W,I}) where {W, I <: Union{Int64,UInt64}} = evfdiv(vone(SVec{W,Float64}), vconvert(SVec{W,Float64}, x))
@inline vinv(x::_Vec{W,I}) where {W, I <: Union{Int32,UInt32}} = evfdiv(vone(_Vec{W,Float32}), vconvert(_Vec{W,Float32}, x))
@inline vinv(x::AbstractStructVec{W,I}) where {W, I <: Union{Int32,UInt32}} = evfdiv(vone(SVec{W,Float32}), vconvert(SVec{W,Float32}, x))
@inline vinv(x::_Vec{W,I}) where {W, I <: Union{Int16,UInt16}} = evfdiv(vone(_Vec{W,Float32}), vconvert(_Vec{W,Float32}, x))
@inline vinv(x::AbstractStructVec{W,I}) where {W, I <: Union{Int16,UInt16}} = evfdiv(vone(SVec{W,Float32}), vconvert(SVec{W,Float32}, x))
@inline Base.inv(v::AbstractStructVec{W,I}) where {W, I<: IntegerTypes} = vinv(v)
@inline Base.sqrt(v::AbstractStructVec{W,I}) where {W, I<: IntegerTypes} = SVec(vsqrt(extract_data(float(v))))

# for accumulating vector results of different sizes.
@generated function vadd(v1::_Vec{_W1,T}, v2::_Vec{_W2,T}) where {_W1,_W2,T}
    W1 = _W1 + 1; W2 = _W2 + 1
    @assert ispow2(W1)
    @assert ispow2(W2)
    W3 = W1
    W4 = W2
    v1e = :v1
    while W3 < W2
        W3 <<= 1
        v1e = Expr(:call, :zeropad, v1e)
    end
    v2e = :v2
    while W4 < W1
        W4 <<= 1
        v2e = Expr(:call, :zeropad, v2e)
    end
    Expr(:block, Expr(:meta, :inline), Expr(:call, :vadd, v1e, v2e))
end
@inline Base.:(+)(x::AbstractStructVec{W1,T}, y::AbstractStructVec{W2,T}) where {W1,W2,T} = SVec(vadd(extract_data(x), extract_data(y)))

@inline Base.abs2(v::SVec) = vmul(v,v)
@inline vabs2(v) = vmul(v,v)
@inline vsum(s::FloatingTypes) = s
@inline vprod(s::FloatingTypes) = s

@inline reduced_add(v::_Vec{W,T}, s::T) where {W,T<:IntegerTypes} = Base.FastMath.add_fast(s, vsum(v))
@inline reduced_add(v::AbstractStructVec{W,T}, s::T) where {W,T<:IntegerTypes} = Base.FastMath.add_fast(s, vsum(v))
@inline reduced_add(s::T, v::_Vec{W,T}) where {W,T<:IntegerTypes} = Base.FastMath.add_fast(s, vsum(v))
@inline reduced_add(s::T, v::AbstractStructVec{W,T}) where {W,T<:IntegerTypes} = Base.FastMath.add_fast(s, vsum(v))

@inline reduced_add(v::AbstractStructVec{W,T}, s::T) where {W,T<:FloatingTypes} = reduced_add(extract_data(v), s)
@inline reduced_add(s::T, v::_Vec{W,T}) where {W,T<:FloatingTypes} = reduced_add(v, s)
@inline reduced_add(s::T, v::AbstractStructVec{W,T}) where {W,T<:FloatingTypes} = reduced_add(extract_data(v), s)

# @inline reduced_add(v1::AbstractSIMDVector{W,T}, v2::AbstractSIMDVector{W,T}) where {W,T} = vadd(v1, v2)
@inline reduced_prod(v::_Vec{W,T}, s::T) where {W,T<:IntegerTypes} = Base.FastMath.mul_fast(s, vprod(v))
@inline reduced_prod(s::T, v::_Vec{W,T}) where {W,T<:IntegerTypes} = Base.FastMath.mul_fast(s, vprod(v))
@inline reduced_prod(v::AbstractStructVec{W,T}, s::T) where {W,T<:IntegerTypes} = Base.FastMath.mul_fast(s, vprod(v))
@inline reduced_prod(s::T, v::AbstractStructVec{W,T}) where {W,T<:IntegerTypes} = Base.FastMath.mul_fast(s, vprod(v))

@inline reduced_prod(s::T, v::_Vec{W,T}) where {W,T<:FloatingTypes} = reduced_prod(v, s)
@inline reduced_prod(v::AbstractStructVec{W,T}, s::T) where {W,T<:FloatingTypes} = reduced_prod(extract_data(v), s)
@inline reduced_prod(s::T, v::AbstractStructVec{W,T}) where {W,T<:FloatingTypes} = reduced_prod(extract_data(v), s)
# @inline reduced_prod(v1::AbstractSIMDVector{W,T}, v2::AbstractSIMDVector{W,T}) where {W,T} = vmul(v1, v2)

@inline reduced_add(v::SVec{W,T}, s::T) where {W,T<:IntegerTypes} = reduced_add(extract_data(v), s)
@inline reduced_prod(v::SVec{W,T}, s::T) where {W,T<:IntegerTypes} = reduced_prod(extract_data(v), s)
@inline reduced_add(v::SVec{W,T}, s::T) where {W,T<:FloatingTypes} = reduced_add(extract_data(v), s)
@inline reduced_prod(v::SVec{W,T}, s::T) where {W,T<:FloatingTypes} = reduced_prod(extract_data(v), s)

@inline reduced_add(v::AbstractSIMDVector{W}, s::T2) where {W,T2} = Base.FastMath.add_fast(s, convert(T2,vsum(v)))
@inline reduced_add(s::T2, v::AbstractSIMDVector{W}) where {W,T2} = Base.FastMath.add_fast(s, convert(T2,vsum(v)))
@inline reduced_add(v1::AbstractSIMDVector{W}, v2::AbstractSIMDVector{W}) where {W} = vadd(v1, v2)
@inline reduced_add(v1::T, v2::T) where {T<:Number} = Base.FastMath.add_fast(v1,v2)
# @inline function reduced_add(v1::T, v2::T) where {T<:Number}
#     @show v1, v2 typeof(v1), typeof(v2)
#     vadd(v1,v2)
# end
@inline reduced_prod(v::AbstractSIMDVector{W}, s::T2) where {W,T2} = Base.FastMath.mul_fast(s, convert(T2,vprod(v)))
@inline reduced_prod(s::T2, v::AbstractSIMDVector{W}) where {W,T2} = Base.FastMath.mul_fast(s, convert(T2,vprod(v)))
@inline reduced_prod(v1::AbstractSIMDVector{W}, v2::AbstractSIMDVector{W}) where {W} = vmul(v1, v2)
@inline reduced_prod(v1::T, v2::T) where {T<:Number} = Base.FastMath.mul_fast(v1,v2)

@inline reduce_to_add(v::_Vec{W,T}, ::T) where {W,T} = vsum(v)
@inline reduce_to_prod(v::_Vec{W,T}, ::T) where {W,T} = vprod(v)
@inline reduce_to_add(v::AbstractStructVec{W,T}, ::T) where {W,T} = vsum(v)
@inline reduce_to_prod(v::AbstractStructVec{W,T}, ::T) where {W,T} = vprod(v)
@inline reduce_to_add(v::AbstractSIMDVector{W}, ::T2) where {W,T2} = convert(T2,vsum(v))
@inline reduce_to_add(v::T, ::T) where {T<:Number} = v
@inline reduce_to_prod(v::AbstractSIMDVector{W}, ::T2) where {W,T2} = convert(T2,vprod(v))
@inline reduce_to_add(v::AbstractSIMDVector{W}, ::AbstractSIMDVector{W}) where {W} = v
@inline reduce_to_prod(v::AbstractSIMDVector{W}, ::AbstractSIMDVector{W}) where {W} = v
@inline reduce_to_prod(s::T, ::T) where {T<:Number} = v


@inline reduced_max(v::_Vec{W,T}, s::T) where {W,T} = max(vmaximum(v), s)
@inline reduced_min(v::_Vec{W,T}, s::T) where {W,T} = min(vminimum(v), s)
@inline reduced_max(v1::_Vec{W,T}, v2::_Vec{W,T}) where {W,T} = vmax(v1, v2)
@inline reduced_min(v1::_Vec{W,T}, v2::_Vec{W,T}) where {W,T} = vmin(v1, v2)
@inline reduce_to_max(v::_Vec{W,T}, s::T) where {W,T} = vmaximum(v)
@inline reduce_to_min(v::_Vec{W,T}, s::T) where {W,T} = vminimum(v)
@inline reduce_to_max(v::_Vec{W,T}, ::_Vec{W,T}) where {W,T} = v
@inline reduce_to_min(v::_Vec{W,T}, ::_Vec{W,T}) where {W,T} = v

@inline reduced_max(v::AbstractStructVec{W,T}, s::T) where {W,T} = max(vmaximum(v), s)
@inline reduced_min(v::AbstractStructVec{W,T}, s::T) where {W,T} = min(vminimum(v), s)
@inline reduced_max(v1::AbstractStructVec{W,T}, v2::AbstractStructVec{W,T}) where {W,T} = vmax(v1, v2)
@inline reduced_min(v1::AbstractStructVec{W,T}, v2::AbstractStructVec{W,T}) where {W,T} = vmin(v1, v2)
@inline reduce_to_max(v::AbstractStructVec{W,T}, s::T) where {W,T} = vmaximum(v)
@inline reduce_to_min(v::AbstractStructVec{W,T}, s::T) where {W,T} = vminimum(v)
@inline reduce_to_max(v::AbstractStructVec{W,T}, ::AbstractStructVec{W,T}) where {W,T} = v
@inline reduce_to_min(v::AbstractStructVec{W,T}, ::AbstractStructVec{W,T}) where {W,T} = v

@inline reduced_max(v::AbstractSIMDVector{W}, s::T2) where {W,T2} = max(convert(T2,vmaximum(v)), s)
@inline reduced_min(v::AbstractSIMDVector{W}, s::T2) where {W,T2} = min(convert(T2,vminimum(v)), s)
@inline reduced_max(v1::AbstractSIMDVector{W}, v2::AbstractStructVec{W,T2}) where {W,T2} = vmax(vconvert(SVec{W,T2},v1), v2)
@inline reduced_min(v1::AbstractSIMDVector{W,}, v2::AbstractStructVec{W,T2}) where {W,T2} = vmin(vconvert(SVec{W,T2},v1), v2)
@inline reduce_to_max(v::AbstractSIMDVector{W}, s::T2) where {W,T2} = convert(T2,vmaximum(v))
@inline reduce_to_min(v::AbstractSIMDVector{W}, s::T2) where {W,T2} = convert(T2,vminimum(v))
@inline reduce_to_max(v::AbstractSIMDVector{W}, ::AbstractStructVec{W,T2}) where {W,T2} = vconvert(SVec{W,T2}, v)
@inline reduce_to_min(v::AbstractSIMDVector{W}, ::AbstractStructVec{W,T2}) where {W,T2} = vconvert(SVec{W,T2}, v)

@inline vnmul(x,y) = vsub(vmul(x, y))
@inline vnsub(x,y) = vsub(y, x)
@inline vmul2(x) = vadd(x,x)
@inline vmul3(x::T) where {T <: Number} = Base.FastMath.mul_fast(T(3), x)
@inline vmul3(v::V) where {W, V <: AbstractSIMDVector{W}} = vmul(vbroadcast(V, T(3)), v)
@inline vadd1(x::T) where {T <: Number} = Base.FastMath.add_fast(x, one(x))
@inline vadd1(v::V) where {W, V <: AbstractSIMDVector{W}} = vadd(vone(V), v)

@generated function addscalar(v::_Vec{_W,T}, s::T) where {_W, T <: Integer}
    W = _W + 1
    typ = "i$(8sizeof(T))"
    vtyp = "<$W x $typ>"
    instrs = String[]
    push!(instrs, "%ie = insertelement $vtyp zeroinitializer, $typ %1, i32 0")
    push!(instrs, "%v = add $vtyp %0, %ie")
    push!(instrs, "ret $vtyp %v")
    quote
        $(Expr(:meta,:inline))
        llvmcall( $(join(instrs,"\n")), NTuple{$W,Core.VecElement{$T}}, Tuple{NTuple{$W,Core.VecElement{$T}},$T}, v, s )
    end
end
@generated function addscalar(v::_Vec{_W,T}, s::T) where {_W, T <: Union{Float16,Float32,Float64}}
    W = _W + 1
    typ = llvmtype(T)
    vtyp = "<$W x $typ>"
    instrs = String[]
    push!(instrs, "%ie = insertelement $vtyp zeroinitializer, $typ %1, i32 0")
    push!(instrs, "%v = fadd $(fastflags(T)) $vtyp %0, %ie")
    push!(instrs, "ret $vtyp %v")
    quote
        $(Expr(:meta,:inline))
        llvmcall( $(join(instrs,"\n")), NTuple{$W,Core.VecElement{$T}}, Tuple{NTuple{$W,Core.VecElement{$T}},$T}, v, s )
    end
end
@inline addscalar(v::SVec, s) = SVec(addscalar(extract_data(v), s))
@inline addscalar(s::T, v::_Vec{W,T}) where {W,T} = addscalar(v, s)
@inline addscalar(s::T, v::AbstractStructVec{W,T}) where {W,T} = addscalar(v, s)
@inline addscalar(a, b) = vadd(a, b)

@generated function mulscalar(v::_Vec{_W,T}, s::T) where {_W, T <: Integer}
    W = _W + 1
    typ = "i$(8sizeof(T))"
    vtyp = "<$W x $typ>"
    instrs = String[]
    push!(instrs, "%ie = insertelement $vtyp $(llvmconst(W, T, 1)), $typ %1, i32 0")
    push!(instrs, "%v = mul $vtyp %0, %ie")
    push!(instrs, "ret $vtyp %v")
    quote
        $(Expr(:meta,:inline))
        llvmcall( $(join(instrs,"\n")), NTuple{$W,Core.VecElement{$T}}, Tuple{NTuple{$W,Core.VecElement{$T}},$T}, v, s )
    end
end
@generated function mulscalar(v::_Vec{_W,T}, s::T) where {_W, T <: Union{Float16,Float32,Float64}}
    W = _W + 1
    typ = llvmtype(T)
    vtyp = "<$W x $typ>"
    instrs = String[]
    push!(instrs, "%ie = insertelement $vtyp $(llvmconst(W, T, 1.0)), $typ %1, i32 0")
    push!(instrs, "%v = fmul $(fastflags(T)) $vtyp %0, %ie")
    push!(instrs, "ret $vtyp %v")
    quote
        $(Expr(:meta,:inline))
        llvmcall( $(join(instrs,"\n")), NTuple{$W,Core.VecElement{$T}}, Tuple{NTuple{$W,Core.VecElement{$T}},$T}, v, s )
    end
end
@inline mulscalar(v::SVec, s) = SVec(mulscalar(extract_data(v), s))
@inline mulscalar(s::T, v::_Vec{W,T}) where {W,T} = mulscalar(v, s)
@inline mulscalar(s::T, v::AbstractStructVec{W,T}) where {W,T} = mulscalar(v, s)
@inline mulscalar(a, b) = vmul(a, b)


function scalar_maxmin(W, ::Type{T}, ismax) where {T}
    comp = if T <: Signed
        LLVM_INS_Int[ismax ? :(>) : :(<)]
    elseif T <: Unsigned
        LLVM_INS_UInt[ismax ? :(>) : :(<)]
    elseif ismax
        "fcmp ogt"
    else
        "fcmp olt"
    end
    basevalue = if T <: Integer
        minmaxzero = ifelse(ismax, typemin(T), typemax(T))
        llvmconst(W, T, minmaxzero)
    else
        llvmconst(W, T, repr(reinterpret(UInt64, ifelse(ismax, -Inf, Inf))))
    end
    typ = llvmtype(T)
    vtyp = "<$W x $typ>"
    instrs = String[]
    push!(instrs, "%ie = insertelement $vtyp $(basevalue), $typ %1, i32 0")
    push!(instrs, "%selection = $comp $vtyp %0, %ie")
    push!(instrs, "%v = select <$W x i1> %selection, $vtyp %0, $vtyp %ie")
    push!(instrs, "ret $vtyp %v")
    instrs
end
@generated function maxscalar(v::_Vec{_W,T}, s::T) where {_W, T}
    W = _W + 1
    instrs = scalar_maxmin(W, T, true)
    quote
        $(Expr(:meta,:inline))
        llvmcall( $(join(instrs,"\n")), NTuple{$W,Core.VecElement{$T}}, Tuple{NTuple{$W,Core.VecElement{$T}},$T}, v, s )
    end
end
@generated function minscalar(v::_Vec{_W,T}, s::T) where {_W, T}
    W = _W + 1
    instrs = scalar_maxmin(W, T, false)
    quote
        $(Expr(:meta,:inline))
        llvmcall( $(join(instrs,"\n")), NTuple{$W,Core.VecElement{$T}}, Tuple{NTuple{$W,Core.VecElement{$T}},$T}, v, s )
    end
end
@inline maxscalar(v::SVec, s) = SVec(maxscalar(extract_data(v), s))
@inline maxscalar(s::T, v::_Vec{W,T}) where {W,T} = maxscalar(v, s)
@inline maxscalar(s::T, v::AbstractStructVec{W,T}) where {W,T} = maxscalar(v, s)
@inline maxscalar(a, b) = vmax(a, b)
@inline minscalar(v::SVec, s) = SVec(minscalar(extract_data(v), s))
@inline minscalar(s::T, v::_Vec{W,T}) where {W,T} = minscalar(v, s)
@inline minscalar(s::T, v::AbstractStructVec{W,T}) where {W,T} = minscalar(v, s)
@inline minscalar(a, b) = vmin(a, b)


@inline Base.literal_pow(::typeof(^), x::AbstractStructVec, ::Val{-2}) = (xi = vinv(x); vmul(xi, xi))
@inline Base.literal_pow(::typeof(^), x::AbstractStructVec, ::Val{-1}) = vinv(x)
@inline Base.literal_pow(::typeof(^), x::V, ::Val{0}) where {V <: AbstractStructVec} = vone(V)
@inline Base.literal_pow(::typeof(^), x::AbstractStructVec, ::Val{1}) = x
@inline Base.literal_pow(::typeof(^), x::AbstractStructVec, ::Val{2}) = vmul(x, x)
@inline Base.literal_pow(::typeof(^), x::AbstractStructVec, ::Val{3}) = vmul(x, vmul(x, x))
@inline function Base.literal_pow(::typeof(^), x::AbstractStructVec, ::Val{4})
    x2 = vmul(x, x)
    vmul(x2, x2)
end
@inline function Base.literal_pow(::typeof(^), x::AbstractStructVec, ::Val{5})
    x2 = vmul(x, x)
    vmul(x, vmul(x2, x2))
end
@inline function Base.literal_pow(::typeof(^), x::AbstractStructVec, ::Val{6})
    x2 = vmul(x, x)
    vmul(x2, vmul(x2, x2))
end
@inline function Base.literal_pow(::typeof(^), x::AbstractStructVec, ::Val{7})
    x2 = vmul(x, x)
    x3 = vmul(x, x2)
    x4 = vmul(x2, x2)
    vmul(x3, x4)
end
@inline function Base.literal_pow(::typeof(^), x::AbstractStructVec, ::Val{8})
    x2 = vmul(x, x)
    x4 = vmul(x2, x2)
    vmul(x4, x4)
end
Base.literal_pow(::typeof(^), x::AbstractStructVec, ::Val{P}) where {P} = x ^ P
# @inline literal_power(x, ::Val{P}) = Base.literal_pow(Base.^, x, Val{P}())

@inline relu(x) = max(x,zero(x))

# These don't work yet.
# @static if Base.libllvm_version ≥ v"10.0.0"
#     @generated function vcolumnwiseload(ptr::Ptr{T}, stride::Integer, ::Val{M}, ::Val{N}) where {T,M,N}
#         typ = llvmtype(T)
#         W = M * N
#         vtyp = "<$W x $typ>"
#         ptyp = "i$(8sizeof(Int))"
#         # suffix = 'r' * M * 'c' * N * typ
#         suffix = "v$(W)$(T <: Integer ? 'i' : 'f')$(8sizeof(T))"
#         instr = "@llvm.matrix.columnwise.load.$suffix"
#         decl = "declare $vtyp $instr($vtyp*, i32, i32, i32)"
#         instrs = String[]
#         push!(instrs, "%ptr = inttoptr $ptyp %0 to $vtyp*")
#         push!(instrs, "%res = call $vtyp $instr($vtyp* %ptr, i32 %1, i32 $M, i32 $N)")
#         push!(instrs, "ret $vtyp %res")
#         quote
#             $(Expr(:meta,:inline))
#             llvmcall(
#                 $((decl, join(instrs, "\n"))),
#                 Vec{$W, $T}, Tuple{Ptr{$T}, UInt32},
#                 ptr, stride % UInt32
#             )
#         end
#     end
#     @generated function vcolumnwisestore!(ptr::Ptr{T}, v::Vec{W,T}, stride::Integer, ::Val{M}, ::Val{N}) where {W,T,M,N}
#         typ = llvmtype(T)
#         W = M * N
#         vtyp = "<$W x $typ>"
#         ptyp = "i$(8sizeof(Int))"
#         # suffix = 'r' * M * 'c' * N * typ
#         suffix = "v$(W)$(T <: Integer ? 'i' : 'f')$(8sizeof(T))"
#         instr = "@llvm.matrix.columnwise.store.$suffix"
#         decl = "declare void $instr($vtyp, $typ*, i32, i32, i32)"
#         instrs = String[]
#         push!(instrs, "%ptr = inttoptr $ptyp %0 to $vtyp*")
#         push!(instrs, "%res = call void $instr($vtyp %1, $vtyp* %ptr, i32 %2, i32 $M, i32 $N)")
#         push!(instrs, "ret void")
#         quote
#             $(Expr(:meta,:inline))
#             llvmcall(
#                 $((decl, join(instrs, "\n"))),
#                 Cvoid, Tuple{Ptr{$T}, Vec{$W,$T}, UInt32},
#                 ptr, v, stride % UInt32
#             )
#         end
#     end
#     @generated function vmatmul(vA::Vec{WA,T}, vB::Vec{WB,T}, ::Val{M}, ::Val{K}, ::Val{N}) where {WA, WB, M, K, N, T}
#         typ = llvmtype(
#         WC = M * N
#         vtypA = "<$WA x $typ>"
#         vtypB = "<$WB x $typ>"
#         vtypC = "<$WC x $typ>"
#         sufftyp = "$(T <: Integer ? 'i' : 'f')$(8sizeof(T))"
#         suffix = "v$(WC)$(sufftyp).v$(WA)$(sufftyp).v$(WB)$(sufftyp)"
#         instr = "@llvm.matrix.multiply.$suffix"
#         decl = "declare $vtypC $instr($vtypA, $vtypB, i32, i32, i32)"
#         instrs = String[]
#         push!(instrs, "%res = call $vtypC $instr($vtypA %0, $vtypB %1, i32 $M, i32 $N, i32 $K)")
#         push!(instrs, "ret $vtypC %res")
#         quote
#             $(Expr(:meta, :inline))
#             llvmcall(
#                 $((decl, join(instrs, "\n"))),
#                 Vec{$WC,$T}, Tuple{Vec{$WA,$T}, Vec{$WB,$T}},
#                 vA, vB
#             )
#         end
#     end
        
    
# # else
# end


