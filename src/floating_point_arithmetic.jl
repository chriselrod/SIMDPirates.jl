
# Floating point arithmetic functions

for op in (
        :+, :-,
        :abs, :ceil, :cos, :exp, :exp2, :floor, :inv, :log, :log10, :log2,
        :round, :sin, :sqrt, :trunc)
    @eval begin
        @inline $op(v1::Vec{N,T}) where {N,T<:FloatingTypes} =
            llvmwrap(Val{$(QuoteNode(op))}, v1)
    end
end
@inline exp10(v1::Vec{N,T}) where {N,T<:FloatingTypes} = broadcast(Vec{N,T}, 10)^v1
@inline sign(v1::Vec{N,T}) where {N,T<:FloatingTypes} =
    vifelse(v1 == broadcast(Vec{N,T}, 0.0), broadcast(Vec{N,T}, 0.0), copysign(broadcast(Vec{N,T}, 1.0), v1))

for op in (:+, :-, :*, :/, :^, :copysign, :max, :min, :rem)
    @eval begin
        @inline $op(v1::Vec{N,T}, v2::Vec{N,T}) where {N,T<:FloatingTypes} =
            llvmwrap(Val{$(QuoteNode(op))}, v1, v2)
    end
end
@inline  ^(v1::Vec{N,T}, x2::Integer) where {N,T<:FloatingTypes} =
    llvmwrap(Val{:powi}, v1, Int(x2))
@inline flipsign(v1::Vec{N,T}, v2::Vec{N,T}) where {N,T<:FloatingTypes} =
    vifelse(signbit(v2), -v1, v1)

for op in (:fma, :muladd)
    @eval begin
        @inline function $op(v1::Vec{N,T},
                v2::Vec{N,T}, v3::Vec{N,T}) where {N,T<:FloatingTypes}
            llvmwrap(Val{$(QuoteNode(op))}, v1, v2, v3)
        end
    end
end

# Type promotion

# Promote scalars of all IntegerTypes to vectors of IntegerTypes, leaving the
# vector type unchanged

for op in (
        :(==), :(!=), :(<), :(<=), :(>), :(>=),
        :&, :|, :⊻, :+, :-, :*, :copysign, :div, :flipsign, :max, :min, :rem)
    @eval begin
        @inline $op(s1::Bool, v2::Vec{N,Bool}) where {N} =
            $op(Vec{N,Bool}(s1), v2)
        @inline $op(s1::IntegerTypes, v2::Vec{N,T}) where {N,T<:IntegerTypes} =
            $op(broadcast(Vec{N,T}, s1), v2)
        @inline $op(v1::Vec{N,T}, s2::IntegerTypes) where {N,T<:IntegerTypes} =
            $op(v1, broadcast(Vec{N,T}, s2))
    end
end
@inline vifelse(c::Vec{N,Bool}, s1::IntegerTypes,
        v2::Vec{N,T}) where {N,T<:IntegerTypes} =
    vifelse(c, broadcast(Vec{N,T}, s1), v2)
@inline vifelse(c::Vec{N,Bool}, v1::Vec{N,T},
        s2::IntegerTypes) where {N,T<:IntegerTypes} =
vifelse(c, v1, broadcast(Vec{N,T}, s2))

for op in (:muladd,)
    @eval begin
        @inline $op(s1::IntegerTypes, v2::Vec{N,T},
                v3::Vec{N,T}) where {N,T<:IntegerTypes} =
            $op(broadcast(Vec{N,T}, s1), v2, v3)
        @inline $op(v1::Vec{N,T}, s2::IntegerTypes,
                v3::Vec{N,T}) where {N,T<:IntegerTypes} =
            $op(v1, broadcast(Vec{N,T}, s2), v3)
        @inline $op(s1::IntegerTypes, s2::IntegerTypes,
                v3::Vec{N,T}) where {N,T<:IntegerTypes} =
            $op(broadcast(Vec{N,T}, s1), broadcast(Vec{N,T}, s2), v3)
        @inline $op(v1::Vec{N,T}, v2::Vec{N,T},
                s3::IntegerTypes) where {N,T<:IntegerTypes} =
            $op(v1, v2, broadcast(Vec{N,T}, s3))
        @inline $op(s1::IntegerTypes, v2::Vec{N,T},
                s3::IntegerTypes) where {N,T<:IntegerTypes} =
            $op(broadcast(Vec{N,T}, s1), v2, broadcast(Vec{N,T}, s3))
        @inline $op(v1::Vec{N,T}, s2::IntegerTypes,
                s3::IntegerTypes) where {N,T<:IntegerTypes} =
            $op(v1, broadcast(Vec{N,T}, s2), broadcast(Vec{N,T}, s3))
    end
end


# Promote scalars of all ScalarTypes to vectors of FloatingTypes, leaving the
# vector type unchanged

for op in (
        :(==), :(!=), :(<), :(<=), :(>), :(>=),
        :+, :-, :*, :/, :^, :copysign, :flipsign, :max, :min, :rem)
    @eval begin
        @inline $op(s1::ScalarTypes, v2::Vec{N,T}) where {N,T<:FloatingTypes} =
            $op(broadcast(Vec{N,T}, s1), v2)
        @inline $op(v1::Vec{N,T}, s2::ScalarTypes) where {N,T<:FloatingTypes} =
            $op(v1, broadcast(Vec{N,T}, s2))
    end
end
@inline vifelse(c::Vec{N,Bool}, s1::ScalarTypes,
        v2::Vec{N,T}) where {N,T<:FloatingTypes} =
    vifelse(c, broadcast(Vec{N,T}, s1), v2)
@inline vifelse(c::Vec{N,Bool}, v1::Vec{N,T},
        s2::ScalarTypes) where {N,T<:FloatingTypes} =
vifelse(c, v1, broadcast(Vec{N,T}, s2))

for op in (:fma, :muladd)
    @eval begin
        @inline $op(s1::ScalarTypes, v2::Vec{N,T},
                v3::Vec{N,T}) where {N,T<:FloatingTypes} =
            $op(broadcast(Vec{N,T}, s1), v2, v3)
        @inline $op(v1::Vec{N,T}, s2::ScalarTypes,
                v3::Vec{N,T}) where {N,T<:FloatingTypes} =
            $op(v1, broadcast(Vec{N,T}, s2), v3)
        @inline $op(s1::ScalarTypes, s2::ScalarTypes,
                v3::Vec{N,T}) where {N,T<:FloatingTypes} =
            $op(broadcast(Vec{N,T}, s1), broadcast(Vec{N,T}, s2), v3)
        @inline $op(v1::Vec{N,T}, v2::Vec{N,T},
                s3::ScalarTypes) where {N,T<:FloatingTypes} =
            $op(v1, v2, broadcast(Vec{N,T}, s3))
        @inline $op(s1::ScalarTypes, v2::Vec{N,T},
                s3::ScalarTypes) where {N,T<:FloatingTypes} =
            $op(broadcast(Vec{N,T}, s1), v2, broadcast(Vec{N,T}, s3))
        @inline $op(v1::Vec{N,T}, s2::ScalarTypes,
                s3::ScalarTypes) where {N,T<:FloatingTypes} =
            $op(v1, broadcast(Vec{N,T}, s2), broadcast(Vec{N,T}, s3))
    end
end

# Reduction operations

# TODO: map, mapreduce

function getneutral(op::Symbol, ::Type{T}) where T
    zs = Dict{Symbol,T}()
    if T <: IntegerTypes
        zs[:&] = ~T(0)
        zs[:|] = T(0)
    end
    zs[:max] = typemin(T)
    zs[:min] = typemax(T)
    zs[:+] = T(0)
    zs[:*] = T(1)
    zs[op]
end

nextpow2(n) = nextpow(2, n)

# We cannot pass in the neutral element via Val{}; if we try, Julia refuses to
# inline this function, which is then disastrous for performance
@generated function llvmwrapreduce(::Type{Val{Op}}, v::Vec{N,T}) where {Op,N,T}
    @assert isa(Op, Symbol)
    z = getneutral(Op, T)
    typ = llvmtype(T)
    decls = []
    instrs = []
    n = N
    nam = "%0"
    nold,n = n,nextpow2(n)
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
        ins = llvmins(Val{Op}, n, T)
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
        Base.llvmcall($((join(decls, "\n"), join(instrs, "\n"))),
            T, Tuple{Vec{N,T}}, v)
    end
end

@inline all(v::Vec{N,T}) where {N,T<:IntegerTypes} = llvmwrapreduce(Val{:&}, v)
@inline any(v::Vec{N,T}) where {N,T<:IntegerTypes} = llvmwrapreduce(Val{:|}, v)
@inline maximum(v::Vec{N,T}) where {N,T<:FloatingTypes} =
    llvmwrapreduce(Val{:max}, v)
@inline minimum(v::Vec{N,T}) where {N,T<:FloatingTypes} =
    llvmwrapreduce(Val{:min}, v)
@inline prod(v::Vec{N,T}) where {N,T} = llvmwrapreduce(Val{:*}, v)
@inline sum(v::Vec{N,T}) where {N,T} = llvmwrapreduce(Val{:+}, v)

@generated function reduce(::Type{Val{Op}}, v::Vec{N,T}) where {Op,N,T}
    @assert isa(Op, Symbol)
    z = getneutral(Op, T)
    stmts = []
    n = N
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

@inline maximum(v::Vec{N,T}) where {N,T<:IntegerTypes} = reduce(Val{:max}, v)
@inline minimum(v::Vec{N,T}) where {N,T<:IntegerTypes} = reduce(Val{:min}, v)
