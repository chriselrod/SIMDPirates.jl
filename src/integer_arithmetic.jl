
# Integer arithmetic functions

for (op,rename) in ((:~,:not), (:+,:add), (:-,:sub))
    @eval begin
        @inline $op(v1::Vec{N,T}) where {N,T<:IntegerTypes} =
            llvmwrap(Val{$(QuoteNode(op))}, v1)
    end
end
# @inline :!(v1::Vec{N,Bool}) where {N} = ~v1
@inline function abs(v1::Vec{N,T}) where {N,T<:IntTypes}
    # s = -broadcast(Vec{N,T}, signbit(v1))
    s = v1 >> Val{8*sizeof(T)}
    # Note: -v1 == ~v1 + 1
    (s ⊻ v1) - s
end
@inline abs(v1::Vec{N,T}) where {N,T<:UIntTypes} = v1
# TODO: Try T(v1>0) - T(v1<0)
#       use a shift for v1<0
#       evaluate v1>0 as -v1<0 ?
@inline sign(v1::Vec{N,T}) where {N,T<:IntTypes} =
    vifelse(v1 == broadcast(Vec{N,T},0), broadcast(Vec{N,T},0),
        vifelse(v1 < broadcast(Vec{N,T},0), broadcast(Vec{N,T},-1), broadcast(Vec{N,T},1)))
@inline sign(v1::Vec{N,T}) where {N,T<:UIntTypes} =
    vifelse(v1 == broadcast(Vec{N,T},0), broadcast(Vec{N,T},0), broadcast(Vec{N,T},1))
@inline signbit(v1::Vec{N,T}) where {N,T<:IntTypes} = v1 < broadcast(Vec{N,T}, 0)
@inline signbit(v1::Vec{N,T}) where {N,T<:UIntTypes} = broadcast(Vec{N,Bool}, false)

for op in (:&, :|, :⊻, :+, :-, :*, :div, :rem)
    @eval begin
        @inline $op(v1::Vec{N,T}, v2::Vec{N,T}) where {N,T<:IntegerTypes} =
            llvmwrap(Val{$(QuoteNode(op))}, v1, v2)
    end
end
@inline copysign(v1::Vec{N,T}, v2::Vec{N,T}) where {N,T<:IntTypes} =
    vifelse(signbit(v2), -abs(v1), abs(v1))
@inline copysign(v1::Vec{N,T}, v2::Vec{N,T}) where {N,T<:UIntTypes} = v1
@inline flipsign(v1::Vec{N,T}, v2::Vec{N,T}) where {N,T<:IntTypes} =
    vifelse(signbit(v2), -v1, v1)
@inline flipsign(v1::Vec{N,T}, v2::Vec{N,T}) where {N,T<:UIntTypes} = v1
@inline max(v1::Vec{N,T}, v2::Vec{N,T}) where {N,T<:IntegerTypes} =
    vifelse(v1>=v2, v1, v2)
@inline min(v1::Vec{N,T}, v2::Vec{N,T}) where {N,T<:IntegerTypes} =
    vifelse(v1>=v2, v2, v1)

@inline function muladd(v1::Vec{N,T}, v2::Vec{N,T},
        v3::Vec{N,T}) where {N,T<:IntegerTypes}
    v1*v2+v3
end
@inline function vfma(v1::Vec{N,T}, v2::Vec{N,T},
        v3::Vec{N,T}) where {N,T<:IntegerTypes}
    v1*v2+v3
end

# TODO: Handle negative shift counts
#       use vifelse
#       ensure vifelse is efficient
for (op,rename) in ((:<<,:left_bitshift), (:>>,:right_bitshift), (:>>>,:uright_bitshift))
    @eval begin
        @inline $rename(v1::Vec{N,T}, ::Type{Val{I}}) where {N,T<:IntegerTypes,I} =
            llvmwrapshift(Val{$(QuoteNode(op))}, v1, Val{I})
        @inline $rename(v1::Vec{N,T}, x2::Unsigned) where {N,T<:IntegerTypes} =
            llvmwrapshift(Val{$(QuoteNode(op))}, v1, x2)
        @inline $rename(v1::Vec{N,T}, x2::Int) where {N,T<:IntegerTypes} =
            llvmwrapshift(Val{$(QuoteNode(op))}, v1, x2)
        @inline $rename(v1::Vec{N,T}, x2::Integer) where {N,T<:IntegerTypes} =
            llvmwrapshift(Val{$(QuoteNode(op))}, v1, x2)
        @inline $rename(v1::Vec{N,T},
                         v2::Vec{N,U}) where {N,T<:IntegerTypes,U<:UIntTypes} =
            llvmwrapshift(Val{$(QuoteNode(op))}, v1, v2)
        @inline $rename(v1::Vec{N,T},
                         v2::Vec{N,U}) where {N,T<:IntegerTypes,U<:IntegerTypes} =
            llvmwrapshift(Val{$(QuoteNode(op))}, v1, v2)
        @inline $rename(x1::T, v2::Vec{N,T}) where {N,T<:IntegerTypes} =
            $rename(broadcast(Vec{N,T}, x1), v2)
    end
end
