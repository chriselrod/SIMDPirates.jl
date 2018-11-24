function _pirate(ex)
    postwalk(ex) do x
        # @show x
        if @capture(x, vadd(vmul(a_, b_), c_)) || @capture(x, vadd(c_, vmul(a_, b_)))
            ea = isa(a, Symbol) ? esc(a) : a
            eb = isa(b, Symbol) ? esc(b) : b
            ec = isa(c, Symbol) ? esc(c) : c
            return :(vmuladd($ea, $eb, $ec))
        elseif @capture(x, vadd(vmul(a_, b_), vmul(c_, d_), e_)) || @capture(x, vadd(vmul(a_, b_), e_, vmul(c_, d_))) || @capture(x, vadd(e_, vmul(a_, b_), vmul(c_, d_)))
            ea = isa(a, Symbol) ? esc(a) : a
            eb = isa(b, Symbol) ? esc(b) : b
            ec = isa(c, Symbol) ? esc(c) : c
            ed = isa(d, Symbol) ? esc(d) : d
            ee = isa(e, Symbol) ? esc(e) : e
            return :(vmuladd($ea, $eb, vmuladd($ec, $ed, $ee)))
        elseif @capture(x, a_ += b_)
            ea = isa(a, Symbol) ? esc(a) : a
            eb = isa(b, Symbol) ? esc(b) : b
            return :($ea = vadd($ea, $eb))
        elseif @capture(x, a_ -= b_)
            ea = isa(a, Symbol) ? esc(a) : a
            eb = isa(b, Symbol) ? esc(b) : b
            return :($ea = vsub($ea, $eb))
        elseif @capture(x, a_ *= b_)
            ea = isa(a, Symbol) ? esc(a) : a
            eb = isa(b, Symbol) ? esc(b) : b
            return :($ea = vmul($ea, $eb))
        elseif @capture(x, a_ /= b_)
            ea = isa(a, Symbol) ? esc(a) : a
            eb = isa(b, Symbol) ? esc(b) : b
            return :($ea = vdiv($ea, $eb))
        elseif isa(x, Symbol) && !occursin("@", string(x))
            return get(VECTOR_SYMBOLS, x, esc(x))
        else
            return x
        end
    end
end


macro pirate(ex) _pirate(ex) end
