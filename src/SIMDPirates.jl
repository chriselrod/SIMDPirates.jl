module SIMDPirates

using VectorizationBase
using VectorizationBase:
    llvmtype, AbstractSIMDVector, SVec, vbroadcast, vzero, vone, _MM, AbstractZeroInitializedPointer,
    AbstractPointer, AbstractInitializedPointer, AbstractStridedPointer, JuliaPointerType
# using MacroTools: prewalk, postwalk
    
export  Vec, SVec, VE,
    @pirate,
    extract_data,
    vbroadcast, vconvert,
    vload,
    vloada,
    vstore!,
    vstorea!,
    shufflevector,
    vifelse,
    vfma, vmuladd,
    vsqrt, rsqrt, vinv, vabs2,
    vadd, vsub, vmul, vfdiv, vsum, vprod,
    vnmul, vnsub,
    vfmadd, vfnmadd, vfmsub, vfnmsub,
    vfmadd_fast, vfnmadd_fast, vfmsub_fast, vfnmsub_fast,
    gather, scatter!,
    addscalar


vecarguments(args) = [isa(arg, Symbol) ? :($arg::Vec{N,T})             : arg for arg ∈ args]
abstractarguments(args) = [isa(arg, Symbol) ? :($arg::AbstractSIMDVector{N,T})    : arg for arg ∈ args]
structvecarguments(args) = [isa(arg, Symbol) ? :($arg::SVec{N,T})     : arg for arg ∈ args]
function vector_args(args)
    vecargs = vecarguments(args)
    abstractargs = abstractarguments(args)
    structargs = structvecarguments(args)
    vecargs, abstractargs, structargs
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
    vecargs, abstractargs, structargs = vector_args(args)
    q = quote
        @inline $rename($(vecargs...)) where {$(R...)} = $body
        @inline $rename($(abstractargs...)) where {$(R...)} = SVec($body)
        @inline $f($(structargs...)) where {$(R...)} = SVec($body)
    end
    esc(q)
end
macro evectordef(rename, expr)
    f, args, R, body = parse_func_decl(expr)
    vecargs, abstractargs, structargs = vector_args(args)
    q = quote
        @inline $rename($(vecargs...)) where {$(R...)} = $body
        @inline $rename($(abstractargs...)) where {$(R...)} = SVec($body)
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
            SVec($(Expr(:call, f, [ Expr(:call, :extract_data, arg) for arg ∈ argnames  ]... )))
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
