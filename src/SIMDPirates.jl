module SIMDPirates

using VectorizationBase
using VectorizationBase:
    llvmtype, AbstractSIMDVector, AbstractStructVec,
    AbstractPointer, AbstractInitializedPointer, AbstractStridedPointer
using MacroTools: prewalk, postwalk
    
export  Vec, SVec, VE,
    @pirate,
    extract_data,
    vbroadcast,
    vload,
    vloada,
    vstore!,
    vstorea!,
    shufflevector,
    vifelse,
    vfma, vmuladd,
    vsqrt, rsqrt, vinv, vabs2,
    vadd, vsub, vmul, vfdiv, vsum, vprod,
    vfmadd, vfnmadd, vfmsub, vfnmsub,
    gather, scatter!


function vector_args(args)
    vecargs =       [isa(arg, Symbol) ? :($arg::Vec{N,T})             : arg for arg ∈ args]
    abstractargs =  [isa(arg, Symbol) ? :($arg::AbstractSIMDVector{N,T})    : arg for arg ∈ args]
    structargs =    [isa(arg, Symbol) ? :($arg::AbstractStructVec{N,T})     : arg for arg ∈ args]
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

include("type_definitions.jl")
include("llvm_utils.jl")
include("llvmwrap.jl")
include("conditionals.jl")
include("integer_arithmetic.jl")
include("floating_point_arithmetic.jl")
include("memory.jl")
include("shufflevector.jl")
include("special.jl")
include("contract_pass.jl")
include("pirate.jl")
include("precompile.jl")
_precompile_()

function __init__()
    _precompile_()
end

end # module
