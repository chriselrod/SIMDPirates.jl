
# Integer arithmetic functions

for op ∈ (:(~), :(+), :(-))
    rename = VECTOR_SYMBOLS[op]
    @eval begin
        @inline $rename(s1::IntegerTypes) = $op($IntegerTypes)
        @inline $rename(v1::Vec{N,T}) where {N,T<:IntegerTypes} =
            llvmwrap(Val{$(QuoteNode(op))}, v1)
    end
end
@inline vnot(s1::Bool) = !s1
@inline vnot(v1::Vec{N,Bool}) where {N} = ~v1
@inline vabs(s1::IntTypes) = abs(s1)
@inline function vabs(v1::Vec{N,T}) where {N,T<:IntTypes}
    # s = -vbroadcast(Vec{N,T}, signbit(v1))
    s = vright_bitshift(v1, Val{8*sizeof(T)}())
    # Note: -v1 == ~v1 + 1
    vsub(vxor(s, v1), s)
end
@inline vabs(v1::Vec{N,T}) where {N,T<:UIntTypes} = v1
# TODO: Try T(v1>0) - T(v1<0)
#       use a shift for v1<0
#       evaluate v1>0 as -v1<0 ?
@inline vsign(s1::IntTypes) = sign(s1)
@inline vsign(v1::Vec{N,T}) where {N,T<:IntTypes} =
    vifelse(vequal(v1, vbroadcast(Vec{N,T},0)), vbroadcast(Vec{N,T},0),
        vifelse(vless(v1, vbroadcast(Vec{N,T},0)), vbroadcast(Vec{N,T},-1), vbroadcast(Vec{N,T},1)))
@inline vsign(v1::Vec{N,T}) where {N,T<:UIntTypes} =
    vifelse(vequal(v1, vbroadcast(Vec{N,T},0)), vbroadcast(Vec{N,T},0), vbroadcast(Vec{N,T},1))
@inline vsignbit(s1::IntegerTypes) = signbit(s1)
@inline vsignbit(v1::Vec{N,T}) where {N,T<:IntTypes} = vless(v1, vbroadcast(Vec{N,T}, 0))
@inline vsignbit(v1::Vec{N,T}) where {N,T<:UIntTypes} = vbroadcast(Vec{N,Bool}, false)

for op ∈ (:(&), :(|), :(⊻), :(+), :(-), :(*), :(÷), :(%) )
    rename = VECTOR_SYMBOLS[op]
    @eval begin
        @inline $rename(s1::IntegerTypes, s2::IntegerTypes) = $op(s1, s2)
        @inline $rename(v1::Vec{N,T}, v2::Vec{N,T}) where {N,T<:IntegerTypes} =
            llvmwrap(Val{$(QuoteNode(op))}, v1, v2)
    end
end

@inline vcopysign(s1::IntegerTypes, s2::IntegerTypes) = copysign(s1, s2)
@inline vcopysign(v1::Vec{N,T}, v2::Vec{N,T}) where {N,T<:IntTypes} =
    vifelse(vsignbit(v2), -abs(v1), abs(v1))
@inline vcopysign(v1::Vec{N,T}, v2::Vec{N,T}) where {N,T<:UIntTypes} = v1

@inline vflipsign(s1::IntegerTypes, s2::IntegerTypes) = flipsign(s1, s2)
@inline vflipsign(v1::Vec{N,T}, v2::Vec{N,T}) where {N,T<:IntTypes} =
    vifelse(vsignbit(v2), -v1, v1)
@inline vflipsign(v1::Vec{N,T}, v2::Vec{N,T}) where {N,T<:UIntTypes} = v1

@inline vmax(s1::IntegerTypes, s2::IntegerTypes) = vmax(s1, s2)
@inline vmax(v1::Vec{N,T}, v2::Vec{N,T}) where {N,T<:IntegerTypes} =
    vifelse(v1>=v2, v1, v2)
@inline vmin(v1::Vec{N,T}, v2::Vec{N,T}) where {N,T<:IntegerTypes} =
    vifelse(v1>=v2, v2, v1)

vmuladd(s1::ScalarTypes, s2::ScalarTypes, s3::ScalarTypes) = muladd(s1,s2,s3)
@inline function vmuladd(v1::Vec{N,T}, v2::Vec{N,T},
        v3::Vec{N,T}) where {N,T<:IntegerTypes}
    vadd(vmul(v1,v2),v3)
end

vfma(s1::ScalarTypes, s2::ScalarTypes, s3::ScalarTypes) = fma(s1,s2,s3)
@inline function vfma(v1::Vec{N,T}, v2::Vec{N,T},
        v3::Vec{N,T}) where {N,T<:IntegerTypes}
    vadd(vmul(v1,v2),v3)
end

# TODO: Handle negative shift counts
#       use vifelse
#       ensure vifelse is efficient
for op ∈ (:(<<), :(>>), :(>>>))
    rename = VECTOR_SYMBOLS[op]
    @eval begin
        @inline $rename(s1::IntegerTypes, s2::IntegerTypes) = $op(s1, s2)
        @inline $rename(v1::Vec{N,T}, ::Val{I}) where {N,T<:IntegerTypes,I} =
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
            $rename(vbroadcast(Vec{N,T}, x1), v2)
    end
end
