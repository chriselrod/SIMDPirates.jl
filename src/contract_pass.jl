
function check_negative(x)
    x isa Expr || return false
    x.head === :call || return false
    length(x.args) == 2 || return false
    a = first(x.args)
    return (a === :(-) || a == :(Base.FastMath.sub_fast))
end

function mulexpr(mulexargs)
    a = (mulexargs[1])::Union{Symbol,Expr}
    b = if length(mulexargs) == 2 # two arg mul
        (mulexargs[2])::Union{Symbol,Expr}
    else
        Expr(:call, :vmul, @view(mulexargs[2:end])...)::Expr
    end
    a, b
end
function append_args_skip!(call, args, i)
    for j ∈ eachindex(args)
        j == i && continue
        push!(call.args, args[j])
    end
    call
end



function recursive_mul_search!(call, argv, cnmul::Bool = false, csub::Bool = false)
    length(argv) < 2 && return length(call.args) == 4, cnmul, csub
    fun = first(argv)
    isadd = fun === :+ || fun === :vadd || fun == :(Base.FastMath.add_fast)
    issub = fun === :- || fun === :vsub || fun == :(Base.FastMath.sub_fast)
    (isadd | issub) || return length(call.args) == 4, cnmul, csub
    exargs = @view(argv[2:end])
    issub && @assert length(exargs) == 2
    for (i,ex) ∈ enumerate(exargs)
        if ex isa Expr && ex.head === :call
            exa = ex.args
            f = first(exa)
            exav = @view(exa[2:end])
            if f === :* || f === :vmul || f == :(Base.FastMath.mul_fast)
                # isnmul = any(check_negative, exav)
                a, b = mulexpr(exav)
                call.args[2] = a
                call.args[3] = b
                if length(exargs) == 2
                    push!(call.args, exargs[3 -  i])
                else
                    push!(call.args, append_args_skip!(Expr(:call, :+), exargs, i))
                end
                if issub
                    csub = i == 1
                    cnmul = !csub
                end
                return true, cnmul, csub
            elseif isadd
                found, cnmul, csub = recursive_mul_search!(call, exa)
                if found
                    if csub
                        call.args[4] = if length(exargs) == 2
                            Expr(:call, :-, exargs[3 - i], call.args[4])
                        else
                            Expr(:call, :-, append_args_skip!(Expr(:call, :+), exargs, i), call.args[4])
                        end
                    else
                        call.args[4] = append_args_skip!(Expr(:call, :+, call.args[4]), exargs, i)
                    end
                    return true, cnmul, false
                end
            elseif issub
                found, cnmul, csub = recursive_mul_search!(call, exa)
                if found
                    if i == 1
                        if csub
                            call.args[4] = Expr(:call, :+, call.args[4], exargs[3 - i])
                        else
                            call.args[4] = Expr(:call, :-, call.args[4], exargs[3 - i])
                        end
                    else
                        cnmul = !cnmul
                        if csub
                            call.args[4] = Expr(:call, :+, exargs[3 - i], call.args[4])
                        else
                            call.args[4] = Expr(:call, :-, exargs[3 - i], call.args[4])
                        end
                        csub = false
                    end
                    return true, cnmul, csub
                end                
            end
        end
    end
    length(call.args) == 4, cnmul, csub
end
                          
function capture_muladd(ex::Expr, mod)
    call = Expr(:call, Symbol(""), Symbol(""), Symbol(""))
    found, nmul, sub = recursive_mul_search!(call, ex.args)
    found || return ex
    if mod === nothing
        call.args[1] = if nmul && sub
            :vfnmsub#_fast
        elseif nmul
            :vfnmadd#_fast
        elseif sub
            :vfmsub#_fast
        else
            :vfmadd#_fast
        end
    else
        call.args[1] = if nmul && sub
            Expr(:(.), mod, QuoteNote(:vfnmsub))#_fast))
        elseif nmul
            Expr(:(.), mod, QuoteNote(:vfnmadd))#_fast))
        elseif sub
            Expr(:(.), mod, QuoteNote(:vfmsub))#_fast))
        else
            Expr(:(.), mod, QuoteNote(:vfmadd))#_fast))
        end
    end
    call
end


contract_pass(x) = x # x will probably be a symbol
function contract_pass(expr::Expr, mod = nothing)::Expr
    prewalk(expr) do ex
        if !(ex isa Expr)
            return ex
        elseif ex.head !== :call
            if ex.head === :(+=)
                call = Expr(:call, :(+))
                append!(call.args, ex.args)
                Expr(:(=), first(ex.args), call)
            elseif ex.head === :(-=)
                call = Expr(:call, :(-))
                append!(call.args, ex.args)
                Expr(:(=), first(ex.args), call)
            elseif ex.head === :(*=)
                call = Expr(:call, :(*))
                append!(call.args, ex.args)
                Expr(:(=), first(ex.args), call)
            elseif ex.head === :(/=)
                call = Expr(:call, :(/))
                append!(call.args, ex.args)
                Expr(:(=), first(ex.args), call)
            else
                ex
            end
        else # ex.head === :call
            return capture_muladd(ex, mod)
        end
    end
end


