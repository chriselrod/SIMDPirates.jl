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
