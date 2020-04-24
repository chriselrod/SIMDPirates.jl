
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
# @inline vexp10(s1::FloatingTypes) = exp10(s1)
# @vectordef vexp10 function Base.exp10(v1) where {W,T <: FloatingTypes}
    # vpow(vbroadcast(Vec{W,T}, 10), extract_data(v1))
# end
# @inline vexp10(v1::AbstractSIMDVector{W,T}) where {W,T<:FloatingTypes} = vpow(vbroadcast(Vec{W,T}, 10), v1)
@inline vsign(s1::FloatingTypes) = sign(s1)
@vectordef vsign function Base.sign(v1) where {W,T<:FloatingTypes}
    vifelse(
        extract_data(v1) == vbroadcast(Vec{W,T}, zero(T)),
        vbroadcast(Vec{W,T}, zero(T)),
        copysign(vbroadcast(Vec{W,T}, one(T)), extract_data(v1))
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
            llvmwrap_notfast(Val{$(QuoteNode(op))}(), extract_data(v1), extract_data(v2))
        end

        @inline function Base.$op(v1::SVec{W,T1}, v2::SVec{W,T2}) where {W,T1,T2}
            T = promote_type(T1, T2)
            $op(vconvert(SVec{W,T}, v1), vconvert(SVec{W,T}, v2))
        end
        @inline function $rename(v1::V1, v2::V2) where {W,T1,T2,V1<:AbstractSIMDVector{W,T1},V2<:AbstractSIMDVector{W,T2}}
            V = promote_vtype(V1, V2)
            $rename(vconvert(V, v1), vconvert(V, v2))
        end
        # @inline function $rename(v1::V1, v2::V2) where {V1,V2}
            # V = vpromote(V1, V2)
            # $rename(vconvert(V, v1), vconvert(V, v2))
        # end
        # @inline $rename(v1::Vec{W,T}, v2::Vec{W,T}) where {W,T<:FloatingTypes} =
        #     llvmwrap(Val{$(QuoteNode(op))}(), v1, v2)
        # @inline $rename(v1::AbstractSIMDVector{W,T}, v2::AbstractSIMDVector{W,T}) where {W,T<:FloatingTypes} =
        #     SVec(llvmwrap(Val{$(QuoteNode(op))}(), extract_data(v1), extract_data(v2)))
        # @inline Base.$op(v1::AbstractSIMDVector{W,T}, v2::AbstractSIMDVector{W,T}) where {W,T<:FloatingTypes} =
        #     SVec(llvmwrap(Val{$(QuoteNode(op))}(), extract_data(v1), extract_data(v2)))
    end
end
@inline vfdiv(v1::Vec{W,Int32}, v2::Vec{W,Int32}) where {W} = vfdiv(vconvert(Vec{W,Float32}, v1), vconvert(Vec{W,Float32}, v2))
@inline vfdiv(v1::Vec{W,Int64}, v2::Vec{W,Int64}) where {W} = vfdiv(vconvert(Vec{W,Float64}, v1), vconvert(Vec{W,Float64}, v2))
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


# let op = :(*)
#     rename = VECTOR_SYMBOLS[op]
#     @eval begin
#         @inline $rename(s1::FloatingTypes, s2::FloatingTypes) = $op(s1, s2)
#         # @inline $rename(v1::Vec{W,T}, v2::Vec{W,T}) where {W,T<:FloatingTypes} =
#         #     VecProduct(v1, v2)
#         @inline $rename(v1::Vec{W,T}, v2::Vec{W,T}) where {W,T<:FloatingTypes} =
#             VecProduct(extract_data(v1), extract_data(v2))
#         @inline $rename(v1::AbstractSIMDVector{W,T}, v2::AbstractSIMDVector{W,T}) where {W,T<:FloatingTypes} =
#             SVecProduct(extract_data(v1), extract_data(v2))
#         @inline Base.$op(v1::SVec{W,T}, v2::SVec{W,T}) where {W,T<:FloatingTypes} =
#             SVecProduct(extract_data(v1), extract_data(v2))
#     end
# end




# @inline vpow(s1::FloatingTypes, x2::Integer) = s1^x2
# @vectordef vpow function Base.:^(v1, x2::Integer) where {W,T<:FloatingTypes}
    # llvmwrap(Val{:powi}, extract_data(v1), Int(x2))
# end
# @inline function vpow(v1::Vec{W,T}, x2::Integer) where {W,T<:FloatingTypes}
#     llvmwrap(Val{:powi}, v1, Int(x2))
# end
# @inline function vpow(v1::SVec{W,T}, x2::Integer) where {W,T<:FloatingTypes}
#     SVec(llvmwrap(Val{:powi}, extract_data(v1), Int(x2)))
# end
# @inline function Base.:^(v1::SVec{W,T}, x2::Integer) where {W,T<:FloatingTypes}
#     SVec(llvmwrap(Val{:powi}, extract_data(v1), Int(x2)))
# end

@inline vflipsign(s1::FloatingTypes, s2::FloatingTypes) = flipsign(s1, s2)

@vectordef vflipsign function Base.flipsign(v1, v2) where {W,T<:FloatingTypes}
    vifelse(vsignbit(extract_data(v2)), vsub(extract_data(v1)), extract_data(v1))
end

@inline vcopysign(v1::Vec{W,T}, v2::Vec{W,U}) where {W,T,U<:Unsigned} = vcopysign(v1, vreinterpret(Vec{W,T}, v2))

# @inline vflipsign(v1::Vec{W,T}, v2::Vec{W,T}) where {W,T<:FloatingTypes} =
#     vifelse(vsignbit(v2), -v1, v1)
# @inline vflipsign(v1::AbstractSIMDVector{W,T}, v2::AbstractSIMDVector{W,T}) where {W,T<:FloatingTypes} =
#     SVec(vifelse(vsignbit(v2), -v1, v1))
# @inline Base.flipsign(v1::SVec{W,T}, v2::SVec{W,T}) where {W,T<:FloatingTypes} =
#     SVec(vifelse(vsignbit(v2), -v1, v1))

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
        # scalar default already set in integer_arithmetic.jl
        # @inline function $rename(v1::Vec{W,T},
        #         v2::Vec{W,T}, v3::Vec{W,T}) where {W,T<:FloatingTypes}
        #     llvmwrap(Val{$(QuoteNode(op))}(), v1, v2, v3)
        # end
        # @inline function $rename(v1::AbstractSIMDVector{W,T},
        #         v2::AbstractSIMDVector{W,T}, v3::AbstractSIMDVector{W,T}) where {W,T<:FloatingTypes}
        #     SVec(llvmwrap(Val{$(QuoteNode(op))}(), extract_data(v1), extract_data(v2), extract_data(v3)))
        # end
        # @inline function Base.$op(v1::SVec{W,T},
        #         v2::SVec{W,T}, v3::SVec{W,T}) where {W,T<:FloatingTypes}
        #     SVec(llvmwrap(Val{$(QuoteNode(op))}(), extract_data(v1), extract_data(v2), extract_data(v3)))
        # end
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
            $rename(vbroadcast(Vec{W,T}, s1), extract_data(v2))
        end
        @vectordef $rename function Base.$op(v1, s2::ScalarTypes) where {W,T<:FloatingTypes}
            $rename(extract_data(v1), vbroadcast(Vec{W,T}, s2))
        end
        
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
        @inline function $rename(s1::T2, v2::Vec{W,T}) where {W,T<:Integer,T2<:FloatingTypes}
            $rename(vbroadcast(Vec{W,T2}, s1), v2)
        end
        @inline function $rename(v1::Vec{W,T}, s2::T2) where {W,T<:Integer,T2<:FloatingTypes}
            $rename(v1, vbroadcast(Vec{W,T2}, s2))
        end
        @vpromote $rename 2
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

# We cannot pass in the neutral element via Val{}; if we try, Julia refuses to
# inline this function, which is then disastrous for performance
@generated function llvmwrapreduce(::Val{Op}, v::Vec{W,T}) where {Op,W,T}
    @assert isa(Op, Symbol)
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
        Base.llvmcall(
            $((join(decls, "\n"), join(instrs, "\n"))),
            $T, Tuple{Vec{$W,$T}}, v
        )
    end
end

@static if Base.libllvm_version >= v"9"
    @generated function vsum(v::Vec{W,T}) where {W,T<:FloatingTypes}
        decls = String[]
        instrs = String[]
        typ = llvmtype(T)
        vtyp = "<$W x $typ>"
        bits = 8sizeof(T)
        ins = "@llvm.experimental.vector.reduce.v2.fadd.f$(bits).v$(W)f$(bits)"
        push!(decls, "declare $(typ) $(ins)($(typ), $(vtyp))")
        push!(instrs, "%res = call fast $typ $ins($typ 0.0, $vtyp %0)")
        push!(instrs, "ret $typ %res")
        quote
            $(Expr(:meta, :inline))
            Base.llvmcall(
                $((join(decls, "\n"), join(instrs, "\n"))),
                $T, Tuple{Vec{$W,$T}}, v
            )
        end
    end
    @generated function vprod(v::Vec{W,T}) where {W,T<:FloatingTypes}
        decls = String[]
        instrs = String[]
        typ = llvmtype(T)
        vtyp = "<$W x $typ>"
        bits = 8sizeof(T)
        ins = "@llvm.experimental.vector.reduce.v2.fmul.f$(bits).v$(W)f$(bits)"
        push!(decls, "declare $(typ) $(ins)($(typ), $(vtyp))")
        push!(instrs, "%res = call fast $typ $ins($typ 1.0, $vtyp %0)")
        push!(instrs, "ret $typ %res")
        quote
            $(Expr(:meta, :inline))
            Base.llvmcall(
                $((join(decls, "\n"), join(instrs, "\n"))),
                $T, Tuple{Vec{$W,$T}}, v
            )
        end
    end
    @generated function reduced_add(v::Vec{W,T}, s::T) where {W,T<:FloatingTypes}
        decls = String[]
        instrs = String[]
        typ = llvmtype(T)
        vtyp = "<$W x $typ>"
        bits = 8sizeof(T)
        ins = "@llvm.experimental.vector.reduce.v2.fadd.f$(bits).v$(W)f$(bits)"
        push!(decls, "declare $(typ) $(ins)($(typ), $(vtyp))")
        push!(instrs, "%res = call fast $typ $ins($typ %1, $vtyp %0)")
        push!(instrs, "ret $typ %res")
        quote
            $(Expr(:meta, :inline))
            Base.llvmcall(
                $((join(decls, "\n"), join(instrs, "\n"))),
                $T, Tuple{Vec{$W,$T},$T}, v, s
            )
        end
    end
    @generated function reduced_prod(v::Vec{W,T}, s::T) where {W,T<:FloatingTypes}
        decls = String[]
        instrs = String[]
        typ = llvmtype(T)
        vtyp = "<$W x $typ>"
        bits = 8sizeof(T)
        ins = "@llvm.experimental.vector.reduce.v2.fmul.f$(bits).v$(W)f$(bits)"
        push!(decls, "declare $(typ) $(ins)($(typ), $(vtyp))")
        push!(instrs, "%res = call fast $typ $ins($typ %1, $vtyp %0)")
        push!(instrs, "ret $typ %res")
        quote
            $(Expr(:meta, :inline))
            Base.llvmcall(
                $((join(decls, "\n"), join(instrs, "\n"))),
                $T, Tuple{Vec{$W,$T},$T}, v, s
            )
        end
    end
    @generated function vsum(v::Vec{W,T}) where {W,T<:IntegerTypes}
        decls = String[]
        instrs = String[]
        typ = llvmtype(T)
        vtyp = "<$W x $typ>"
        bits = 8sizeof(T)
        ins = "@llvm.experimental.vector.reduce.add.v$(W)i$(bits)"
        push!(decls, "declare $(typ) $(ins)($(vtyp))")
        push!(instrs, "%res = call $typ $ins($vtyp %0)")
        push!(instrs, "ret $typ %res")
        quote
            $(Expr(:meta, :inline))
            Base.llvmcall(
                $((join(decls, "\n"), join(instrs, "\n"))),
                $T, Tuple{Vec{$W,$T}}, v
            )
        end
    end
    @generated function vprod(v::Vec{W,T}) where {W,T<:IntegerTypes}
        decls = String[]
        instrs = String[]
        typ = llvmtype(T)
        vtyp = "<$W x $typ>"
        bits = 8sizeof(T)
        ins = "@llvm.experimental.vector.reduce.mul.v$(W)i$(bits)"
        push!(decls, "declare $(typ) $(ins)($(vtyp))")
        push!(instrs, "%res = call $typ $ins($vtyp %0)")
        push!(instrs, "ret $typ %res")
        quote
            $(Expr(:meta, :inline))
            Base.llvmcall(
                $((join(decls, "\n"), join(instrs, "\n"))),
                $T, Tuple{Vec{$W,$T}}, v
            )
        end
    end
    @generated function vsub(v::Vec{W,T}) where {W,T<:FloatingTypes}
        typ = llvmtype(T)
        vtyp = "<$W x $typ>"
        instrs = "%res = fneg fast $vtyp %0\nret $vtyp %res"
        quote
            $(Expr(:meta, :inline))
            Base.llvmcall( $instrs, Vec{$W,$T}, Tuple{Vec{$W,$T}}, v )
        end
    end
    Base.:(-)(v::SVec{W,T}) where {W,T} = SVec(vsub(extract_data(v)))
    vsub(v::SVec{W,T}) where {W,T} = SVec(vsub(extract_data(v)))
else
    # @generated function vsub(v::Vec{W,T}) where {W,T<:FloatingTypes}
    #     typ = llvmtype(T)
    #     vtyp = "<$W x $typ>"
    #     instrs = "%res = fsub fast $vtyp zeroinitializer, %0\nret $vtyp %res"
    #     quote
    #         $(Expr(:meta, :inline))
    #         Base.llvmcall( $instrs, Vec{$W,$T}, Tuple{Vec{$W,$T}}, v )
    #     end
    # end
    # Base.:(-)(v::SVec{W,T}) where {W,T} = SVec(vsub(extract_data(v)))
    # vsub(v::SVec{W,T}) where {W,T} = SVec(vsub(extract_data(v)))
    rename = VECTOR_SYMBOLS[:-]
    @eval begin
        @vectordef $rename function Base.:(-)(v1) where {W,T<:FloatingTypes}
            llvmwrap(Val{:(-)}(), extract_data(v1))
        end
    end
end
vsub(x::FloatingTypes) = Base.FastMath.sub_fast(x)

for (name, rename, op) ∈ ((:(Base.all),:vall,:&), (:(Base.any),:vany,:|),
                                    (:(Base.maximum), :vmaximum, :max), (:(Base.minimum), :vminimum, :min),
                                    (:(Base.sum),:vsum,:+), (:(Base.prod),:vprod,:*))
    @eval begin
        @inline $rename(v::AbstractSIMDVector{W,T}) where {W,T} = llvmwrapreduce(Val{$(QuoteNode(op))}(), extract_data(v))
        @inline $name(v::SVec{W,T}) where {W,T} = llvmwrapreduce(Val{$(QuoteNode(op))}(), extract_data(v))
    end
end

# @inline vall(v::Vec{W,T}) where {W,T<:IntegerTypes} = llvmwrapreduce(Val{:&}, v)
# @inline vany(v::Vec{W,T}) where {W,T<:IntegerTypes} = llvmwrapreduce(Val{:|}, v)
# @inline vmaximum(v::Vec{W,T}) where {W,T<:FloatingTypes} =
#     llvmwrapreduce(Val{:max}, v)
# @inline vminimum(v::Vec{W,T}) where {W,T<:FloatingTypes} =
#     llvmwrapreduce(Val{:min}, v)
# @inline vprod(v::Vec{W,T}) where {W,T} = llvmwrapreduce(Val{:*}, v)
# @inline vsum(v::Vec{W,T}) where {W,T} = llvmwrapreduce(Val{:+}, v)

@generated function vreduce(::Val{Op}, v::Vec{W,T}) where {Op,W,T}
    @assert isa(Op, Symbol)
    z = getneutral(Op, T)
    stmts = String[]
    n = W
    push!(stmts, :($(Symbol(:v,n)) = v))
    nold,n = n,nextpow2(n)
    if n > nold
        push!(stmts,
            :($(Symbol(:v,n)) = $(Expr(:tuple,
                [:($(Symbol(:v,nold))[$i]) for i in 1:nold]...,
                [:(VE($z)) for i in nold+1:n]...))))
    end
    while n > 1
        nold,n = n, div(n, 2)
        push!(stmts,
            :($(Symbol(:v,n,"lo")) = $(Expr(:tuple,
                [:($(Symbol(:v,nold))[$i]) for i in 1:n]...,))))
        push!(stmts,
            :($(Symbol(:v,n,"hi")) = $(Expr(:tuple,
                [:($(Symbol(:v,nold))[$i]) for i in n+1:nold]...))))
        push!(stmts,
            :($(Symbol(:v,n)) =
                $Op($(Symbol(:v,n,"lo")), $(Symbol(:v,n,"hi")))))
    end
    push!(stmts, :(v1[1].value))
    Expr(:block, Expr(:meta, :inline), stmts...)
end

@vectordef vmaximum function Base.maximum(v) where {W,T<:IntegerTypes}
    vreduce(Val{:max}, extract_data(v))
end
@vectordef vminimum function Base.minimum(v) where {W,T<:IntegerTypes}
    vreduce(Val{:min}, extract_data(v))
end

# @inline vmaximum(v::Vec{W,T}) where {W,T<:IntegerTypes} = vreduce(Val{:max}, v)
# @inline vminimum(v::Vec{W,T}) where {W,T<:IntegerTypes} = vreduce(Val{:min}, v)

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

@generated function vfmadd(v1::Vec{W,T}, v2::Vec{W,T}, v3::Vec{W,T}) where {W,T<:FloatingTypes}
    typ = llvmtype(T)
    vtyp = "<$W x $typ>"
    # ins = "@llvm.fma.v$(W)f$(8*sizeof(T))"
    ins = "@llvm.fmuladd.v$(W)f$(8*sizeof(T))"
    decls = "declare $vtyp $ins($vtyp, $vtyp, $vtyp)"
    instrs = "%res = call $vtyp $ins($vtyp %0, $vtyp %1, $vtyp %2)\nret $vtyp %res"
    quote
        $(Expr(:meta, :inline))
        Base.llvmcall(
            $((decls,instrs)),
            Vec{$W,$T}, Tuple{Vec{$W,$T}, Vec{$W,$T}, Vec{$W,$T}},
            v1, v2, v3
        )
    end
end
vfmadd(a::Number, b::Number, c::Number) = muladd(a, b, c)
vfnmadd(a::Number, b::Number, c::Number) = muladd(-a, b, c)
vfmsub(a::Number, b::Number, c::Number) = muladd(a, b, -c)
vfnmsub(a::Number, b::Number, c::Number) = -muladd(a, b, c)
@inline vfmadd(a::Vec{W,T}, b::Vec{W,T}, c::Vec{W,T}) where {W,T <: Integer} = vmuladd(a, b, c)
@inline vfnmadd(a::Vec{W,T}, b::Vec{W,T}, c::Vec{W,T}) where {W,T} = vfmadd(vsub(a), b, c)
@inline vfmsub(a::Vec{W,T}, b::Vec{W,T}, c::Vec{W,T}) where {W,T} = vfmadd(a, b, vsub(c))
@inline vfnmsub(a::Vec{W,T}, b::Vec{W,T}, c::Vec{W,T}) where {W,T} = vsub(vfmadd(a, b, c))
@vpromote vfmadd 3
@vpromote vfnmadd 3
@vpromote vfmsub 3
@vpromote vfnmsub 3

@generated function vfmadd_fast(v1::Vec{W,T}, v2::Vec{W,T}, v3::Vec{W,T}) where {W,T<:FloatingTypes}
    typ = llvmtype(T)
    vtyp = "<$W x $typ>"
    # ins = "@llvm.fmuladd.v$(W)f$(8*sizeof(T))"
    # decls = "declare $vtyp $ins($vtyp, $vtyp, $vtyp)"
    # instrs = "%res = call fast $vtyp $ins($vtyp %0, $vtyp %1, $vtyp %2)\nret $vtyp %res"
    # I don't really want reassoc on the fmul part, so this is what I'm going with.
    instrs = """
    %prod = fmul nnan ninf nsz arcp contract $vtyp %0, %1
    %res = fadd fast $vtyp %prod, %2
    ret $vtyp %res
    """
    quote
        $(Expr(:meta, :inline))
        Base.llvmcall(
            $instrs, Vec{$W,$T}, Tuple{Vec{$W,$T}, Vec{$W,$T}, Vec{$W,$T}}, v1, v2, v3
        )
    end
end
vfmadd_fast(a::Number, b::Number, c::Number) = Base.FastMath.add_fast(Base.FastMath.mul_fast(a, b), c)
vfnmadd_fast(a::Number, b::Number, c::Number) = Base.FastMath.sub_fast(c, Base.FastMath.mul_fast(a, b))
vfmsub_fast(a::Number, b::Number, c::Number) = Base.FastMath.sub_fast(Base.FastMath.mul_fast(a, b), c)
vfnmsub_fast(a::Number, b::Number, c::Number) = Base.FastMath.sub_fast(Base.FastMath.add_fast(Base.FastMath.mul_fast(a, b), c))
@inline vfmadd_fast(a::Vec{W,T}, b::Vec{W,T}, c::Vec{W,T}) where {W,T <: Integer} = vmuladd(a, b, c)
@inline vfnmadd_fast(a::Vec{W,T}, b::Vec{W,T}, c::Vec{W,T}) where {W,T} = vfmadd_fast(vsub(a), b, c)
@inline vfmsub_fast(a::Vec{W,T}, b::Vec{W,T}, c::Vec{W,T}) where {W,T} = vfmadd_fast(a, b, vsub(c))
@inline vfnmsub_fast(a::Vec{W,T}, b::Vec{W,T}, c::Vec{W,T}) where {W,T} = vsub(vfmadd_fast(a, b, c))
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

# Lowers to same split mul-add llvm as the 
# @inline vfmadd(a, b, c) = vadd(vmul( a, b), c)
# definition, so I wont bother implementing
# vfnmadd, vfmsub, and vfnmsub
# in this manner.

# @inline Base.:*(a::IntegerTypes, b::SVec{W,T}) where {W,T} = SVec{W,T}(a) * b
# @inline Base.:*(a::T, b::SVec{W,<:IntegerTypes}) where {W,T<:FloatingTypes} = SVec(vmul(vbroadcast(Vec{W,T}, a), vconvert(Vec{W,T}, extract_data(b))))

# const Vec{W,T} = NTuple{W,Core.VecElement{T}}
@inline vfmadd231(a::Vec, b, c) = vfmadd(a, b, c)
@inline vfnmadd231(a::Vec, b, c) = vfnmadd(a, b, c)
@inline vfmsub231(a::Vec, b, c) = vfmsub(a, b, c)
@inline vfnmsub231(a::Vec, b, c) = vfnmsub(a, b, c)
@inline vfmadd231(a::Number, b::Number, c::Number) = Base.FastMath.add_fast(Base.FastMath.mul_fast(a, b), c)
@inline vfnmadd231(a::Number, b::Number, c::Number) = Base.FastMath.sub_fast(c, Base.FastMath.mul_fast(a, b))
@inline vfmsub231(a::Number, b::Number, c::Number) = Base.FastMath.sub_fast(Base.FastMath.mul_fast(a, b), c)
@inline vfnmsub231(a::Number, b::Number, c::Number) = Base.FastMath.sub_fast(Base.FastMath.sub_fast(c), Base.FastMath.mul_fast(a, b))
@inline function vifelse(f::typeof(vfmadd231), m::AbstractMask{W}, a, b, c) where {W}
    vifelse(m, vfmadd_fast(a, b, c), c)
end
@inline function vifelse(f::typeof(vfnmadd231), m::AbstractMask{W}, a, b, c) where {W}
    vifelse(m, vfnmadd_fast(a, b, c), c)
end
@inline function vifelse(f::typeof(vfmsub231), m::AbstractMask{W}, a, b, c) where {W}
    vifelse(m, vfmsub_fast(a, b, c), c)
end
@inline function vifelse(f::typeof(vfnmsub231), m::AbstractMask{W}, a, b, c) where {W}
    vifelse(m, vfnmsub_fast(a, b, c), c)
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
                    Base.llvmcall($vfmadd_str, Vec{$W,$T}, Tuple{Vec{$W,$T},Vec{$W,$T},Vec{$W,$T}}, a, b, c)
                end
                @inline function vfnmadd231(a::Vec{$W,$T}, b::Vec{$W,$T}, c::Vec{$W,$T})
                    Base.llvmcall($vfnmadd_str, Vec{$W,$T}, Tuple{Vec{$W,$T},Vec{$W,$T},Vec{$W,$T}}, a, b, c)
                end
                @inline function vfmsub231(a::Vec{$W,$T}, b::Vec{$W,$T}, c::Vec{$W,$T})
                    Base.llvmcall($vfmsub_str, Vec{$W,$T}, Tuple{Vec{$W,$T},Vec{$W,$T},Vec{$W,$T}}, a, b, c)
                end
                @inline function vfnmsub231(a::Vec{$W,$T}, b::Vec{$W,$T}, c::Vec{$W,$T})
                    Base.llvmcall($vfnmsub_str, Vec{$W,$T}, Tuple{Vec{$W,$T},Vec{$W,$T},Vec{$W,$T}}, a, b, c)
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
#             Base.llvmcall(
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
@inline vinv(x::IntegerTypes) = vinv(float(x))
@inline vinv(x::Vec{W,I}) where {W, I <: Union{Int64,UInt64}} = evfdiv(vone(Vec{W,Float64}), vconvert(Vec{W,Float64}, x))
@inline vinv(x::AbstractSIMDVector{W,I}) where {W, I <: Union{Int64,UInt64}} = evfdiv(vone(SVec{W,Float64}), vconvert(SVec{W,Float64}, x))
@inline vinv(x::Vec{W,I}) where {W, I <: Union{Int32,UInt32}} = evfdiv(vone(Vec{W,Float32}), vconvert(Vec{W,Float32}, x))
@inline vinv(x::AbstractSIMDVector{W,I}) where {W, I <: Union{Int32,UInt32}} = evfdiv(vone(SVec{W,Float32}), vconvert(SVec{W,Float32}, x))
@inline vinv(x::Vec{W,I}) where {W, I <: Union{Int16,UInt16}} = evfdiv(vone(Vec{W,Float32}), vconvert(Vec{W,Float32}, x))
@inline vinv(x::AbstractSIMDVector{W,I}) where {W, I <: Union{Int16,UInt16}} = evfdiv(vone(SVec{W,Float32}), vconvert(SVec{W,Float32}, x))
@inline Base.inv(v::AbstractStructVec{W,I}) where {W, I<: IntegerTypes} = vinv(v)
@inline Base.sqrt(v::AbstractStructVec{W,I}) where {W, I<: IntegerTypes} = vsqrt(float(v))

# for accumulating vector results of different sizes.
@generated function vadd(v1::Vec{W1,T}, v2::Vec{W2,T}) where {W1,W2,T}
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

@inline reduced_add(v::AbstractSIMDVector{W,T}, s::T) where {W,T} = Base.FastMath.add_fast(s, vsum(v))
@inline reduced_add(s::T, v::AbstractSIMDVector{W,T}) where {W,T} = Base.FastMath.add_fast(s, vsum(v))
# @inline reduced_add(v1::AbstractSIMDVector{W,T}, v2::AbstractSIMDVector{W,T}) where {W,T} = vadd(v1, v2)
@inline reduced_prod(v::AbstractSIMDVector{W,T}, s::T) where {W,T} = Base.FastMath.mul_fast(s, vprod(v))
@inline reduced_prod(s::T, v::AbstractSIMDVector{W,T}) where {W,T} = Base.FastMath.mul_fast(s, vprod(v))
# @inline reduced_prod(v1::AbstractSIMDVector{W,T}, v2::AbstractSIMDVector{W,T}) where {W,T} = vmul(v1, v2)

@inline reduced_add(v::SVec{W,T}, s::T) where {W,T} = reduced_add(extract_data(v), s)
@inline reduced_prod(v::SVec{W,T}, s::T) where {W,T} = reduced_prod(extract_data(v), s)
@inline reduced_add(v::AbstractSIMDVector{W,T1}, s::T2) where {W,T1,T2} = Base.FastMath.add_fast(s, convert(T2,vsum(v)))
@inline reduced_add(s::T2, v::AbstractSIMDVector{W,T1}) where {W,T1,T2} = Base.FastMath.add_fast(s, convert(T2,vsum(v)))
@inline reduced_add(v1::AbstractSIMDVector{W,T1}, v2::AbstractSIMDVector{W,T2}) where {W,T1,T2} = vadd(v1, v2)
@inline reduced_add(v1::T, v2::T) where {T<:Number} = Base.FastMath.add_fast(v1,v2)
@inline reduced_prod(v::AbstractSIMDVector{W,T1}, s::T2) where {W,T1,T2} = Base.FastMath.mul_fast(s, convert(T2,vprod(v)))
@inline reduced_prod(s::T2, v::AbstractSIMDVector{W,T1}) where {W,T1,T2} = Base.FastMath.mul_fast(s, convert(T2,vprod(v)))
@inline reduced_prod(v1::AbstractSIMDVector{W,T1}, v2::AbstractSIMDVector{W,T2}) where {W,T1,T2} = vmul(v1, v2)
@inline reduced_prod(v1::T, v2::T) where {T<:Number} = Base.FastMath.mul_fast(v1,v2)

@inline reduce_to_add(v::AbstractSIMDVector{W,T}, ::T) where {W,T} = vsum(v)
@inline reduce_to_prod(v::AbstractSIMDVector{W,T}, ::T) where {W,T} = vprod(v)
@inline reduce_to_add(v::AbstractSIMDVector{W,T1}, ::T2) where {W,T1,T2} = convert(T2,vsum(v))
@inline reduce_to_add(v::T, ::T) where {T<:Number} = v
@inline reduce_to_prod(v::AbstractSIMDVector{W,T1}, ::T2) where {W,T1,T2} = convert(T2,vprod(v))
@inline reduce_to_add(v::AbstractSIMDVector{W,T}, ::AbstractSIMDVector{W}) where {W,T} = v
@inline reduce_to_prod(v::AbstractSIMDVector{W,T}, ::AbstractSIMDVector{W}) where {W,T} = v
@inline reduce_to_prod(s::T, ::T) where {T<:Number} = v


# @inline reduced_all(u1::Unsigned, u1::Unsigned) where {W,T} = vall(v) & s
# @inline reduced_any(u::Unsigned, b::Bool) where {W,T} = b || vany(v)
# @inline reduced_all(v::AbstractSIMDVector, s::T) where {W,T} = vall(v) & s
# @inline reduced_any(v::AbstractSIMDVector, s::T) where {W,T} = s || vany(v)
@inline reduced_max(v::AbstractSIMDVector{W,T}, s::T) where {W,T} = max(vmaximum(v), s)
@inline reduced_min(v::AbstractSIMDVector{W,T}, s::T) where {W,T} = min(vminimum(v), s)
@inline reduced_max(v1::AbstractSIMDVector{W,T}, v2::AbstractSIMDVector{W,T}) where {W,T} = vmax(v1, v2)
@inline reduced_min(v1::AbstractSIMDVector{W,T}, v2::AbstractSIMDVector{W,T}) where {W,T} = vmin(v1, v2)
@inline reduce_to_max(v::AbstractSIMDVector{W,T}, s::T) where {W,T} = vmaximum(v)
@inline reduce_to_min(v::AbstractSIMDVector{W,T}, s::T) where {W,T} = vminimum(v)
@inline reduce_to_max(v::AbstractSIMDVector{W,T}, ::AbstractSIMDVector{W,T}) where {W,T} = v
@inline reduce_to_min(v::AbstractSIMDVector{W,T}, ::AbstractSIMDVector{W,T}) where {W,T} = v

@inline reduced_max(v::AbstractSIMDVector{W,T1}, s::T2) where {W,T1,T2} = max(convert(T2,vmaximum(v)), s)
@inline reduced_min(v::AbstractSIMDVector{W,T1}, s::T2) where {W,T1,T2} = min(convert(T2,vminimum(v)), s)
@inline reduced_max(v1::AbstractSIMDVector{W,T1}, v2::AbstractSIMDVector{W,T2}) where {W,T1,T2} = vmax(vconvert(Vec{W,T2},v1), v2)
@inline reduced_min(v1::AbstractSIMDVector{W,T1}, v2::AbstractSIMDVector{W,T2}) where {W,T1,T2} = vmin(vconvert(Vec{W,T2},v1), v2)
@inline reduce_to_max(v::AbstractSIMDVector{W,T1}, s::T2) where {W,T1,T2} = convert(T2,vmaximum(v))
@inline reduce_to_min(v::AbstractSIMDVector{W,T1}, s::T2) where {W,T1,T2} = convert(T2,vminimum(v))
@inline reduce_to_max(v::AbstractSIMDVector{W,T1}, ::AbstractSIMDVector{W,T2}) where {W,T1,T2} = vconvert(Vec{W,T2}, v)
@inline reduce_to_min(v::AbstractSIMDVector{W,T1}, ::AbstractSIMDVector{W,T2}) where {W,T1,T2} = vconvert(Vec{W,T2}, v)

@inline vnmul(x,y) = vsub(vmul(x, y))
@inline vnsub(x,y) = vsub(y, x)
@inline vmul2(x) = vadd(x,x)
@inline vmul3(x::T) where {T <: Number} = Base.FastMath.mul_fast(T(3), x)
@inline vmul3(v::V) where {W, T, V <: AbstractSIMDVector{W,T}} = vmul(vbroadcast(V, T(3)), v)
@inline vadd1(x::T) where {T <: Number} = Base.FastMath.add_fast(x, one(x))
@inline vadd1(v::V) where {W, T, V <: AbstractSIMDVector{W,T}} = vadd(vbroadcast(V, one(T)), v)

@generated function addscalar(v::Vec{W,T}, s::T) where {W, T <: Integer}
    typ = "i$(8sizeof(T))"
    vtyp = "<$W x $typ>"
    instrs = String[]
    push!(instrs, "%ie = insertelement $vtyp zeroinitializer, $typ %1, i32 0")
    push!(instrs, "%v = add $vtyp %0, %ie")
    push!(instrs, "ret $vtyp %v")
    quote
        $(Expr(:meta,:inline))
        Base.llvmcall( $(join(instrs,"\n")), NTuple{$W,Core.VecElement{$T}}, Tuple{NTuple{$W,Core.VecElement{$T}},$T}, v, s )
    end
end
@generated function addscalar(v::Vec{W,T}, s::T) where {W, T <: Union{Float16,Float32,Float64}}
    typ = llvmtype(T)
    vtyp = "<$W x $typ>"
    instrs = String[]
    push!(instrs, "%ie = insertelement $vtyp zeroinitializer, $typ %1, i32 0")
    push!(instrs, "%v = fadd fast $vtyp %0, %ie")
    push!(instrs, "ret $vtyp %v")
    quote
        $(Expr(:meta,:inline))
        Base.llvmcall( $(join(instrs,"\n")), NTuple{$W,Core.VecElement{$T}}, Tuple{NTuple{$W,Core.VecElement{$T}},$T}, v, s )
    end
end
@inline addscalar(v::SVec, s) = addscalar(extract_data(v), s)
@inline addscalar(s::T, v::Union{Vec{W,T},SVec{W,T}}) where {W,T} = addscalar(v, s)
@inline addscalar(a, b) = vadd(a, b)

@generated function mulscalar(v::Vec{W,T}, s::T) where {W, T <: Integer}
    typ = "i$(8sizeof(T))"
    vtyp = "<$W x $typ>"
    instrs = String[]
    push!(instrs, "%ie = insertelement $vtyp $(llvmconst(W, T, 1)), $typ %1, i32 0")
    push!(instrs, "%v = mul $vtyp %0, %ie")
    push!(instrs, "ret $vtyp %v")
    quote
        $(Expr(:meta,:inline))
        Base.llvmcall( $(join(instrs,"\n")), NTuple{$W,Core.VecElement{$T}}, Tuple{NTuple{$W,Core.VecElement{$T}},$T}, v, s )
    end
end
@generated function mulscalar(v::Vec{W,T}, s::T) where {W, T <: Union{Float16,Float32,Float64}}
    typ = llvmtype(T)
    vtyp = "<$W x $typ>"
    instrs = String[]
    push!(instrs, "%ie = insertelement $vtyp $(llvmconst(W, T, 1.0)), $typ %1, i32 0")
    push!(instrs, "%v = fmul fast $vtyp %0, %ie")
    push!(instrs, "ret $vtyp %v")
    quote
        $(Expr(:meta,:inline))
        Base.llvmcall( $(join(instrs,"\n")), NTuple{$W,Core.VecElement{$T}}, Tuple{NTuple{$W,Core.VecElement{$T}},$T}, v, s )
    end
end
@inline mulscalar(v::SVec, s) = mulscalar(extract_data(v), s)
@inline mulscalar(s::T, v::Union{Vec{W,T},SVec{W,T}}) where {W,T} = mulscalar(v, s)
@inline mulscalar(a, b) = vmul(a, b)

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

