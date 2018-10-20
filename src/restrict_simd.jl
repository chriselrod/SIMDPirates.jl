using CpuId
using MacroTools: prewalk
const REGISTER_SIZE = simdbytes()

mask_expr(W, r) = :($(Expr(:tuple, [i > r ? Core.VecElement{Bool}(false) : Core.VecElement{Bool}(true) for i ∈ 1:W]...)))

function vectorize_body(N, T, n, body)
    T_size = sizeof(T)
    W = REGISTER_SIZE ÷ T_size
    while W > 2N
        W >>= 1
    end
    WT = W * T_size
    Q, r = divrem(N, W) #Assuming Mₖ is a multiple of W
    QQ, Qr = divrem(Q, 4)
    if r > 0
        Qr += 1
        Q += 1
    end
    # unroll the remainder iteration
    # so that whenever Q >= 4, we will always have at least
    # 4 operations scheduled at a time.
    if QQ > 0 && Qr > 0 && Qr < 4 # if r > 0, Qr may equal 4
        QQ -= 1
        Qr += 4
    end
    V = Vec{W,T}


    # body = _pirate(body)

    # indexed_expressions = Dict{Symbol,Expr}()
    indexed_expressions = Dict{Symbol,Symbol}()
    # itersym = esc(gensym(:iter))
    # itersym = esc(:iter)
    itersym = :iter
    # walk the expression, searching for all get index patterns.
    # these will be replaced with
    # Plan: definition of q will create pointers
    main_body = quote end
    for b ∈ body
        push!(main_body.args, _pirate(prewalk(b) do x
            if @capture(x, A_[i_] = B_)
                if A ∉ keys(indexed_expressions)
                    # pA = esc(gensym(A))
                    # pA = esc(Symbol(:p,A))
                    pA = Symbol(:p,A)
                    indexed_expressions[A] = pA
                else
                    pA = indexed_expressions[A]
                end
                eB = isa(B, Symbol) ? esc(B) : B
                if i == n
                    return :(vstore($eB, $pA + $itersym))
                else
                    ei = isa(i, Symbol) ? esc(i) : i
                    eA = isa(A, Symbol) ? esc(A) : A
                    return :(vstore($eB, $pA + $ei*sizeof(eltype($eA))))
                end
            elseif @capture(x, A_[i_])
                if A ∉ keys(indexed_expressions)
                    # pA = esc(gensym(A))
                    # pA = esc(Symbol(:p,A))
                    pA = Symbol(:p,A)
                    indexed_expressions[A] = pA
                else
                    pA = indexed_expressions[A]
                end
                if i == n
                    return :(vload($V, $pA + $itersym))
                else
                    # when loading something not indexed by the loop variable,
                    # we assume that the intension is to broadcast
                    ei = isa(i, Symbol) ? esc(i) : i
                    return :(vbroadcast($V, unsafe_load($pA + $ei*sizeof(eltype($A)))))
                end
            else
                return x
            end
        end))# |> x -> (@show(x), _pirate(x)))
    end

    q = quote end
    for (sym, psym) ∈ indexed_expressions
        push!(q.args, :( $(esc(psym)) = $(esc(Base.pointer))($(esc(sym))) ))
    end

    if QQ > 0
        push!(q.args,
        quote
            for i ∈ 0:$(4W*sizeof(T)):$((QQ-1)*4W*sizeof(T))
                $itersym = i
                $main_body
                $itersym = i + $(W*sizeof(T))
                $main_body
                $itersym = i + $(2W*sizeof(T))
                $main_body
                $itersym = i + $(3W*sizeof(T))
                $main_body
            end
        end)
    end
    for i ∈ 0:Qr-1
        push!(q.args,
        quote
            $itersym = $(QQ*4W*sizeof(T) + i*4W*sizeof(T))
            $main_body
        end)
    end
    if r > 0
        mask = mask_expr(W, r)
        iter = Q * W * sizeof(T)
        r_body = quote end
        for b ∈ body
            push!(r_body.args, _pirate(prewalk(b) do x
                if @capture(x, A_[i_] = B_)
                    if A ∉ keys(indexed_expressions)
                        # pA = esc(gensym(A))
                        # pA = esc(Symbol(:p,A))
                        pA = Symbol(:p,A)
                        indexed_expressions[A] = pA
                    else
                        pA = indexed_expressions[A]
                    end
                    eB = isa(B, Symbol) ? esc(B) : B
                    if i == n
                        return :(vstore($eB, $pA + $iter, $mask))
                    else
                        ei = isa(i, Symbol) ? esc(i) : i
                        eA = isa(A, Symbol) ? esc(A) : A
                        return :(vstore($eB, $pA + $ei*sizeof(eltype($eA)), $mask))
                    end
                elseif @capture(x, A_[i_])
                    if A ∉ keys(indexed_expressions)
                        # pA = esc(gensym(A))
                        # pA = esc(Symbol(:p,A))
                        pA = Symbol(:p,A)
                        indexed_expressions[A] = pA
                    else
                        pA = indexed_expressions[A]
                    end
                    if i == n
                        return :(vload($V, $pA + $iter, $mask))
                    else
                        # when loading something not indexed by the loop variable,
                        # we assume that the intension is to broadcast
                        ei = isa(i, Symbol) ? esc(i) : i
                        return :(vbroadcast($V, unsafe_load($pA + $ei*sizeof(eltype($A)), $mask)))
                    end
                else
                    return x
                end
            end))
        end
        push!(q.args, r_body)
    end
    q
end

macro restrict_simd(expr)
    if @capture(expr, for n_ ∈ 1:N_ body__ end)
        q = vectorize_body(N, Float64, n, body)
    # elseif @capture(expr, for n_ ∈ 1:N_ body__ end)
    #     q = vectorize_body(N, element_type(body)
    else
        throw("Could not match loop expression.")
    end
    q
end
macro restrict_simd(type, expr)
    if @capture(expr, for n_ ∈ 1:N_ body__ end)
        q = vectorize_body(N, type, n, body)
    else
        throw("Could not match loop expression.")
    end
    q
end
