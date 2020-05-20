# Vector shuffles

function shufflevector_instrs(N, T, I, two_operands)
    typ = llvmtype(T)
    vtyp2 = vtyp1 = "<$N x $typ>"
    M = length(I)
    vtyp3 = "<$M x i32>"
    vtypr = "<$M x $typ>"
    mask = "<" * join(map(x->string("i32 ", x), I), ", ") * ">"
    instrs = String[]
    v2 = two_operands ? "%1" : "undef"
    push!(instrs, "%res = shufflevector $vtyp1 %0, $vtyp2 $v2, $vtyp3 $mask")
    push!(instrs, "ret $vtypr %res")
    return M, String[], instrs
end


@generated function shufflevector(
    v1::_Vec{_W,T}, v2::_Vec{_W,T}, ::Val{I}
) where {_W,T,I}
    W = _W + 1
    M, decls, instrs = shufflevector_instrs(W, T, I, true)
    quote
        $(Expr(:meta, :inline))
        Vec{$M,T}(Base.llvmcall($((join(decls, "\n"), join(instrs, "\n"))),
            Vec{$M,T},
            Tuple{Vec{W,T}, Vec{W,T}},
            v1, v2))
    end
end
@inline function shufflevector(v1::SVec{W,T}, v2::SVec{W,T}, ::Val{I}) where {W,T,I}
    SVec(shufflevector(extract_data(v1), extract_data(v2), Val(I)))
end

@generated function shufflevector(v1::_Vec{_W,T}, ::Val{I}) where {_W,T,I}
    W = _W + 1
    M, decls, instrs = shufflevector_instrs(W, T, I, false)
    quote
        $(Expr(:meta, :inline))
        Base.llvmcall(
            $((join(decls, "\n"), join(instrs, "\n"))),
            Vec{$M,T}, Tuple{Vec{$W,T}}, v1
        )
    end
end
@inline function shufflevector(v1::SVec{W,T}, ::Val{I}) where {W,T,I}
    SVec(shufflevector(extract_data(v1), Val(I)))
end

@inline rotate_vector_left(v::AbstractSIMDVector{W}) where {W} = shufflevector(v, Val(ntuple(i -> (i % W), Val(W))))

