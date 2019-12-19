function horner(x, pu...)
    t = gensym(:t)
    ex = p[end]
    for i ∈ length(p)-1:-1:1
        ex = :(SIMDPirates.vmuladd($t, $ex, $(p[i])))
    end
    Expr(:block, :($t = $x), ex)
end

function _pirate(ex)
    postwalk(contract_pass(ex, :SIMDPirates)) do x
        if x isa Symbol
            f = get(VECTOR_SYMBOLS, x, x)
            if f === x
                return x
            else
                return Expr(:(.), :SIMDPirates, QuoteNode(f))
            end
        end
        x isa Expr || return x
        xexpr::Expr = x
        xexpr.head === :call || return x
        f = first(xexpr.args)
        if f == :(Base.FastMath.add_fast)
            vf = :vadd
        elseif f == :(Base.FastMath.sub_fast)
            vf = :vsub
        elseif f == :(Base.FastMath.mul_fast)
            vf = :vmul
        elseif f == :(Base.FastMath.div_fast)
            vf = :vfdiv
        elseif f == :(Base.FastMath.sqrt)
            vf = :vsqrt
        elseif f == :(Base.Math.muladd)
            vf = :vmuladd
        else
            return xexpr
        end
    end |> esc
end


macro pirate(ex) _pirate(ex) end


struct vBitArray
    ptr::Ptr{UInt64}
end
@inline VectorizationBase.vectorizable(b::BitArray) = vBitArray(pointer(b.chunks))
@inline function Base.:+(i, b::vBitArray)
    vBitArray((i >> 3) + b.ptr)
end
@inline function Base.:+(b::vBitArray, i)
    vBitArray((i >> 3) + b.ptr)
end
@inline function Base.:-(b::vBitArray, i)
    vBitArray(b.ptr - (i >> 3))
end
"""
The method on BitArrays is a hack, ignoring the type passed to Vec, to better support the LoopVectorization.jl implementation.
"""
@generated function vload(::Type{Vec{N,T}}, b::vBitArray, i::Integer) where {N,T}
    N < 8 && throw("Bit array vectors with $N < 8 not yet supported.")
    utype = mask_type(N)
    quote
        $(Expr(:meta, :inline))
        unsafe_load(Base.unsafe_convert(Ptr{$utype}, b.ptr + (i >> 3)))
    end
end
@generated function vload(::Type{Vec{N,T}}, b::vBitArray) where {N,T}
    N < 8 && throw("Bit array vectors with $N < 8 not yet supported.")
    utype = mask_type(N)
    quote
        $(Expr(:meta, :inline))
        unsafe_load(Base.unsafe_convert(Ptr{$utype}, b.ptr))
    end
end
@generated function vload(::Type{SVec{N,T}}, b::vBitArray) where {N,T}
    N < 8 && throw("Bit array vectors with $N < 8 not yet supported.")
    utype = mask_type(N)
    quote
        $(Expr(:meta, :inline))
        unsafe_load(Base.unsafe_convert(Ptr{$utype}, b.ptr))
    end
end
@generated function vload(::Type{SVec{N,T}}, b::vBitArray, i) where {N,T}
    N < 8 && throw("Bit array vectors with $N < 8 not yet supported.")
    utype = mask_type(N)
    quote
        $(Expr(:meta, :inline))
        unsafe_load(Base.unsafe_convert(Ptr{$utype}, b.ptr + (i >> 3)))
    end
end

"""
Masks on vectorizable bit arrays are currently ignored.
"""
@generated function vload(::Type{Vec{N,T}}, b::vBitArray, i::Integer, mask::Union{<:Unsigned,Vec{N,Bool}}) where {N,T}
    N < 8 && throw("Bit array vectors with $N < 8 not yet supported.")
    utype = mask_type(N)
    quote
        $(Expr(:meta, :inline))
        unsafe_load(Base.unsafe_convert(Ptr{$utype}, b.ptr + (i >> 3)))
    end
end
@generated function vload(::Type{Vec{N,T}}, b::vBitArray, mask::Union{<:Unsigned,Vec{N,Bool}}) where {N,T}
    N < 8 && throw("Bit array vectors with $N < 8 not yet supported.")
    utype = mask_type(N)
    quote
        $(Expr(:meta, :inline))
        unsafe_load(Base.unsafe_convert(Ptr{$utype}, b.ptr))
    end
end
@generated function vload(::Type{SVec{N,T}}, b::vBitArray, mask::Union{<:Unsigned,Vec{N,Bool}}) where {N,T}
    N < 8 && throw("Bit array vectors with $N < 8 not yet supported.")
    utype = mask_type(N)
    quote
        $(Expr(:meta, :inline))
        unsafe_load(Base.unsafe_convert(Ptr{$utype}, b.ptr))
    end
end
@generated function vload(::Type{SVec{N,T}}, b::vBitArray, i, mask::Union{<:Unsigned,Vec{N,Bool}}) where {N,T}
    N < 8 && throw("Bit array vectors with $N < 8 not yet supported.")
    utype = mask_type(N)
    quote
        $(Expr(:meta, :inline))
        unsafe_load(Base.unsafe_convert(Ptr{$utype}, b.ptr + (i >> 3)))
    end
end
