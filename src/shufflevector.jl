# Vector shuffles

function shufflevector_instrs(N, T, I, two_operands)
    typ = llvmtype(T)
    vtyp2 = vtyp1 = "<$N x $typ>"
    M = length(I)
    vtyp3 = "<$M x i32>"
    vtypr = "<$M x $typ>"
    mask = "<" * join(map(x->string("i32 ", x), I), ", ") * ">"
    instrs = []
    v2 = two_operands ? "%1" : "undef"
    push!(instrs, "%res = shufflevector $vtyp1 %0, $vtyp2 $v2, $vtyp3 $mask")
    push!(instrs, "ret $vtypr %res")
    return M, [], instrs
end


@generated function shufflevector(v1::Vec{N,T}, v2::Vec{N,T},
                                  ::Val{I}) where {N,T,I}
    M, decls, instrs = shufflevector_instrs(N, T, I, true)
    quote
        $(Expr(:meta, :inline))
        Vec{$M,T}(Base.llvmcall($((join(decls, "\n"), join(instrs, "\n"))),
            Vec{$M,T},
            Tuple{Vec{N,T}, Vec{N,T}},
            v1, v2))
    end
end
@inline function shufflevector(v1::AbstractStructVec{N,T}, v2::AbstractStructVec{N,T}, ::Val{I}) where {N,T,I}
    SVec(shufflevector(extract_data(v1), extract_data(v2), Val(I)))
end

@generated function shufflevector(v1::Vec{N,T}, ::Val{I}) where {N,T,I}
    M, decls, instrs = shufflevector_instrs(N, T, I, false)
    quote
        $(Expr(:meta, :inline))
        Base.llvmcall($((join(decls, "\n"), join(instrs, "\n"))),
            Vec{$M,T},
            Tuple{Vec{$N,T}},
            v1)
    end
end
@inline function shufflevector(v1::AbstractStructVec{N,T}, ::Val{I}) where {N,T,I}
    SVec(shufflevector(extract_data(v1), Val(I)))
end
