function horner(x, pu...)
    t = gensym(:t)
    ex = p[end]
    for i ∈ length(p)-1:-1:1
        ex = :(SIMDPirates.vmuladd($t, $ex, $(p[i])))
    end
    Expr(:block, :($t = $x), ex)
end

spfunc(f) = Expr(:(.), :SIMDPirates, QuoteNode(f))
function ret_pirate(ex::Expr, x::Symbol, i::Int)
    f = get(VECTOR_SYMBOLS, x, x)
    f === x || (ex.args[i] = spfunc(f))
    nothing
end
ret_pirate!(::Expr, ::Any, ::Int) = nothing
function ret_pirate!(ex::Expr, x::Expr, i::Int)
    length(x.args) == 0 && return
    f = first(x.args)
    x.args[1] = if f == :(Base.FastMath.add_fast)
        spfunc(:vadd)
    elseif f == :(Base.FastMath.sub_fast)
        spfunc(:vsub)
    elseif f == :(Base.FastMath.mul_fast)
        spfunc(:vmul)
    elseif f == :(Base.FastMath.div_fast)
        spfunc(:vfdiv)
    elseif f == :(Base.FastMath.sqrt)
        spfunc(:vsqrt)
    elseif f == :(Base.Math.muladd)
        spfunc(:vmuladd)
    else
        f
    end
    _pirate!(expr, mod)
    nothing
end
    

function _pirate!(expr, mod)
    for (i,ex) ∈ enumerate(expr.args)
        ret_pirate!(expr, ex, i)
    end
    esc(expr)
end

macro pirate(ex) _pirate!(contract_pass!(ex, Expr(:(.), mod, QuoteNode(:SIMDPirates))), Symbol(__module__)) end


# struct vBitArray
#     ptr::Ptr{UInt64}
# end
# @inline VectorizationBase.vectorizable(b::BitArray) = vBitArray(pointer(b.chunks))
# @inline function Base.:+(i, b::vBitArray)
#     vBitArray((i >> 3) + b.ptr)
# end
# @inline function Base.:+(b::vBitArray, i)
#     vBitArray((i >> 3) + b.ptr)
# end
# @inline function Base.:-(b::vBitArray, i)
#     vBitArray(b.ptr - (i >> 3))
# end
# """
# The method on BitArrays is a hack, ignoring the type passed to Vec, to better support the LoopVectorization.jl implementation.
# """
# @generated function vload(::Type{Vec{N,T}}, b::vBitArray, i::Integer) where {N,T}
#     N < 8 && throw("Bit array vectors with $N < 8 not yet supported.")
#     utype = mask_type(N)
#     quote
#         $(Expr(:meta, :inline))
#         unsafe_load(Base.unsafe_convert(Ptr{$utype}, b.ptr + (i >> 3)))
#     end
# end
# @generated function vload(::Type{Vec{N,T}}, b::vBitArray) where {N,T}
#     N < 8 && throw("Bit array vectors with $N < 8 not yet supported.")
#     utype = mask_type(N)
#     quote
#         $(Expr(:meta, :inline))
#         unsafe_load(Base.unsafe_convert(Ptr{$utype}, b.ptr))
#     end
# end
# @generated function vload(::Type{SVec{N,T}}, b::vBitArray) where {N,T}
#     N < 8 && throw("Bit array vectors with $N < 8 not yet supported.")
#     utype = mask_type(N)
#     quote
#         $(Expr(:meta, :inline))
#         unsafe_load(Base.unsafe_convert(Ptr{$utype}, b.ptr))
#     end
# end
# @generated function vload(::Type{SVec{N,T}}, b::vBitArray, i) where {N,T}
#     N < 8 && throw("Bit array vectors with $N < 8 not yet supported.")
#     utype = mask_type(N)
#     quote
#         $(Expr(:meta, :inline))
#         unsafe_load(Base.unsafe_convert(Ptr{$utype}, b.ptr + (i >> 3)))
#     end
# end

# """
# Masks on vectorizable bit arrays are currently ignored.
# """
# @generated function vload(::Type{Vec{N,T}}, b::vBitArray, i::Integer, mask::Union{<:Unsigned,Vec{N,Bool}}) where {N,T}
#     N < 8 && throw("Bit array vectors with $N < 8 not yet supported.")
#     utype = mask_type(N)
#     quote
#         $(Expr(:meta, :inline))
#         unsafe_load(Base.unsafe_convert(Ptr{$utype}, b.ptr + (i >> 3)))
#     end
# end
# @generated function vload(::Type{Vec{N,T}}, b::vBitArray, mask::Union{<:Unsigned,Vec{N,Bool}}) where {N,T}
#     N < 8 && throw("Bit array vectors with $N < 8 not yet supported.")
#     utype = mask_type(N)
#     quote
#         $(Expr(:meta, :inline))
#         unsafe_load(Base.unsafe_convert(Ptr{$utype}, b.ptr))
#     end
# end
# @generated function vload(::Type{SVec{N,T}}, b::vBitArray, mask::Union{<:Unsigned,Vec{N,Bool}}) where {N,T}
#     N < 8 && throw("Bit array vectors with $N < 8 not yet supported.")
#     utype = mask_type(N)
#     quote
#         $(Expr(:meta, :inline))
#         unsafe_load(Base.unsafe_convert(Ptr{$utype}, b.ptr))
#     end
# end
# @generated function vload(::Type{SVec{N,T}}, b::vBitArray, i, mask::Union{<:Unsigned,Vec{N,Bool}}) where {N,T}
#     N < 8 && throw("Bit array vectors with $N < 8 not yet supported.")
#     utype = mask_type(N)
#     quote
#         $(Expr(:meta, :inline))
#         unsafe_load(Base.unsafe_convert(Ptr{$utype}, b.ptr + (i >> 3)))
#     end
# end
