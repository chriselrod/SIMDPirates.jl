function horner(x, p...)
    t = gensym(:t)
    ex = p[end]
    for i ∈ length(p)-1:-1:1
        ex = :(SIMDPirates.vmuladd($t, $ex, $(p[i])))
    end
    Expr(:block, :($t = $x), ex)
end

function _pirate(ex)
    postwalk(ex) do x
        # @show x
        # if @capture(x, SIMDPirates.vadd(SIMDPirates.vmul(a_, b_), c_)) || @capture(x, SIMDPirates.vadd(c_, SIMDPirates.vmul(a_, b_)))
        #     return :(SIMDPirates.vmuladd($a, $b, $c))
        # elseif @capture(x, SIMDPirates.vadd(SIMDPirates.vmul(a_, b_), SIMDPirates.vmul(c_, d_), e_)) || @capture(x, SIMDPirates.vadd(SIMDPirates.vmul(a_, b_), e_, SIMDPirates.vmul(c_, d_))) || @capture(x, SIMDPirates.vadd(e_, SIMDPirates.vmul(a_, b_), SIMDPirates.vmul(c_, d_)))
        #     return :(SIMDPirates.vmuladd($a, $b, SIMDPirates.vmuladd($c, $d, $e)))
        # elseif @capture(x, a_ += b_)
        if @capture(x, a_ += b_)
            return :($a = SIMDPirates.vadd($a, $b))
        elseif @capture(x, a_ -= b_)
            return :($a = SIMDPirates.vsub($a, $b))
        elseif @capture(x, a_ *= b_)
            return :($a = SIMDPirates.vmul($a, $b))
        elseif @capture(x, a_ /= b_)
            return :($a = SIMDPirates.vdiv($a, $b))
        elseif @capture(x, @horner a__)
            return horner(a...)
        elseif @capture(x, Base.Math.muladd(a_, b_, c_))
            return :( SIMDPirates.vmuladd($a, $b, $c) )
        elseif isa(x, Symbol) && !occursin("@", string(x))
            if x ∈ keys(VECTOR_SYMBOLS)
                return :(SIMDPirates.$(VECTOR_SYMBOLS[x]))
            else
                return x
            end
            # return get(VECTOR_SYMBOLS, x, x)
            # return :(extract_data($(get(VECTOR_SYMBOLS, x, x))))
        else
            return x
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
function unsigned_type(N)
    if N <= 8
        return :UInt8
    elseif N <= 16
        return :UInt16
    elseif N <= 32
        return :UInt32
    elseif N <= 64
        return :UInt64
    elseif N <= 128
        return :UInt128
    else
        throw("Vectors of length $N > 128 not yet supported.")
    end
end
"""
The method on BitArrays is a hack, ignoring the type passed to Vec, to better support the LoopVectorization.jl implementation.
"""
@generated function vload(::Type{Vec{N,T}}, b::vBitArray, i::Integer) where {N,T}
    N < 8 && throw("Bit array vectors with $N < 8 not yet supported.")
    utype = unsigned_type(N)
    quote
        $(Expr(:meta, :inline))
        unsafe_load(Base.unsafe_convert(Ptr{$utype}, b.ptr + (i >> 3)))
    end
end
@generated function vload(::Type{Vec{N,T}}, b::vBitArray) where {N,T}
    N < 8 && throw("Bit array vectors with $N < 8 not yet supported.")
    utype = unsigned_type(N)
    quote
        $(Expr(:meta, :inline))
        unsafe_load(Base.unsafe_convert(Ptr{$utype}, b.ptr))
    end
end
@generated function vload(::Type{SVec{N,T}}, b::vBitArray) where {N,T}
    N < 8 && throw("Bit array vectors with $N < 8 not yet supported.")
    utype = unsigned_type(N)
    quote
        $(Expr(:meta, :inline))
        unsafe_load(Base.unsafe_convert(Ptr{$utype}, b.ptr))
    end
end
@generated function vload(::Type{SVec{N,T}}, b::vBitArray, i) where {N,T}
    N < 8 && throw("Bit array vectors with $N < 8 not yet supported.")
    utype = unsigned_type(N)
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
    utype = unsigned_type(N)
    quote
        $(Expr(:meta, :inline))
        unsafe_load(Base.unsafe_convert(Ptr{$utype}, b.ptr + (i >> 3)))
    end
end
@generated function vload(::Type{Vec{N,T}}, b::vBitArray, mask::Union{<:Unsigned,Vec{N,Bool}}) where {N,T}
    N < 8 && throw("Bit array vectors with $N < 8 not yet supported.")
    utype = unsigned_type(N)
    quote
        $(Expr(:meta, :inline))
        unsafe_load(Base.unsafe_convert(Ptr{$utype}, b.ptr))
    end
end
@generated function vload(::Type{SVec{N,T}}, b::vBitArray, mask::Union{<:Unsigned,Vec{N,Bool}}) where {N,T}
    N < 8 && throw("Bit array vectors with $N < 8 not yet supported.")
    utype = unsigned_type(N)
    quote
        $(Expr(:meta, :inline))
        unsafe_load(Base.unsafe_convert(Ptr{$utype}, b.ptr))
    end
end
@generated function vload(::Type{SVec{N,T}}, b::vBitArray, i, mask::Union{<:Unsigned,Vec{N,Bool}}) where {N,T}
    N < 8 && throw("Bit array vectors with $N < 8 not yet supported.")
    utype = unsigned_type(N)
    quote
        $(Expr(:meta, :inline))
        unsafe_load(Base.unsafe_convert(Ptr{$utype}, b.ptr + (i >> 3)))
    end
end
