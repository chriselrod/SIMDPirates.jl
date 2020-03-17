module SIMDPirates

using VectorizationBase
using VectorizationBase:
    llvmtype, AbstractSIMDVector, SVec, vbroadcast, vzero, vone, _MM,
    AbstractPointer, AbstractInitializedPointer, AbstractStridedPointer, JuliaPointerType, AbstractStructVec
import VectorizationBase: vload, vstore!, vnoaliasstore!


export  Vec, SVec, VE, _MM, stridedpointer,
    extract_data, vbroadcast, vconvert,
    vload, vloada, vstore!, vstorea!, vstorent!, shufflevector,
    vifelse, Mask,
    vfma, vmuladd,
    vsqrt, rsqrt, vinv, vabs2,
    vadd, vsub, vmul, vfdiv, vsum, vprod,
    vnmul, vnsub,
    vfmadd, vfnmadd, vfmsub, vfnmsub,
    vfmadd_fast, vfnmadd_fast, vfmsub_fast, vfnmsub_fast,
    gather, scatter!,
    addscalar,
    @pirate


vecarguments(args) = [isa(arg, Symbol) ? :($arg::Vec{W,T})             : arg for arg ∈ args]
abstractarguments(args) = [isa(arg, Symbol) ? :($arg::AbstractSIMDVector{W,T})    : arg for arg ∈ args]
function structvecarguments(args)
    asa = Expr[]
    sva = Expr[]
    rna = Expr[]
    for arg ∈ args
        if isa(arg, Symbol)
            rarg = gensym(arg)
            push!(asa, Expr(:(::), rarg, :(AbstractStructVec{W,T})))
            push!(sva, Expr(:(::), rarg, :(SVec{W,T})))
            push!(rna, Expr(:(=), arg, Expr(:call, :extract_data, rarg)))
        elseif length(arg.args) == 2
            arg = copy(arg)
            narg = arg.args[1]
            rarg = gensym(narg)
            arg.args[1] = rarg
            push!(asa, arg)
            push!(sva, arg)
            push!(rna, Expr(:(=), narg, Expr(:call, :vconvert, :(Vec{W,T}), rarg)))
        else
            push!(asa, arg)
            push!(sva, arg)
        end
    end
    asa, sva, rna
end
function vector_args(args)
    vecargs = vecarguments(args)
    abstractargs = abstractarguments(args)
    asargs, structargs, renamedargs = structvecarguments(args)
    vecargs, abstractargs, asargs, structargs, renamedargs
end

function parse_func_decl(expr::Expr)
    @assert (expr.head === :(=) || expr.head === :function) "Failed to match function declaration: $expr."
    w = expr.args[1]
    body = expr.args[2]
    WT = @view(w.args[2:end])
    func = first(w.args)
    f = first(func.args)
    args = @view(func.args[2:end])
    f, args, WT, body
end
macro vectordef(rename, expr)
    f, args, R, body = parse_func_decl(expr)
    vecargs, abstractargs, asargs, structargs, renamedargs = vector_args(args)
    push!(renamedargs, Expr(:call, Expr(:curly, :SVec, :W), body))
    q = quote
        @inline $rename($(vecargs...)) where {$(R...)} = $body
        @inline $rename($(abstractargs...)) where {$(R...)} = SVec{W}($body)
        @inline function $f($(asargs...)) where {$(R...)}
            $(renamedargs...)
        end
        @inline function $f($(structargs...)) where {$(R...)}
            $(renamedargs...)
        end
    end
    esc(q)
end
macro evectordef(rename, expr)
    f, args, R, body = parse_func_decl(expr)
    vecargs, abstractargs, structargs = vector_args(args)
    q = quote
        @inline $rename($(vecargs...)) where {$(R...)} = $body
        @inline $rename($(abstractargs...)) where {$(R...)} = SVec{W}($body)
    end
    esc(q)
end

macro vpromote(f, N)
    argnames = [Symbol(:v_, n) for n ∈ 1:N]
    types = [Symbol(:V, n) for n ∈ 1:N]
    args = [Expr(:(::), a, t) for (a,t) ∈ zip(argnames, types)]
    svecargs = [Expr(:(::), a, Expr(:curly, :SVec, :W, :T)) for a ∈ argnames]
    q = quote
        @inline function $f($(args...)) where {$(types...)}
            V = promote_vtype($(types...))
            $(Expr(:call, f, [ Expr(:call, :vconvert, :V, arg) for arg ∈ argnames  ]... ))
        end
        @inline function $f($(svecargs...)) where {W,T}
            SVec{W}($(Expr(:call, f, [ Expr(:call, :extract_data, arg) for arg ∈ argnames  ]... )))
        end
    end
    esc(q)
end

include("type_definitions.jl")
include("llvm_utils.jl")
include("llvmwrap.jl")
include("conditionals.jl")
include("integer_arithmetic.jl")
include("floating_point_arithmetic.jl")
include("vrange.jl")
include("memory.jl")
include("shufflevector.jl")
include("special.jl")
include("contract_pass.jl")
include("pirate.jl")
include("arithmeticwithconsts.jl")
include("precompile.jl")
_precompile_()


end # module
