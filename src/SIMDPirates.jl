module SIMDPirates

using VectorizationBase
using VectorizationBase:
    llvmtype, AbstractSIMDVector, AbstractStructVec, vbroadcast, vzero, vone,
    AbstractPointer, AbstractInitializedPointer, AbstractStridedPointer, JuliaPointerType
using MacroTools: prewalk, postwalk
    
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
    gather, scatter!


vecarguments(args) = [isa(arg, Symbol) ? :($arg::Vec{N,T})             : arg for arg ∈ args]
abstractarguments(args) = [isa(arg, Symbol) ? :($arg::AbstractSIMDVector{N,T})    : arg for arg ∈ args]
structvecarguments(args) = [isa(arg, Symbol) ? :($arg::AbstractStructVec{N,T})     : arg for arg ∈ args]
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

macro vpromote(expr)
    f = first(expr.args)::Symbol
    args = convert(Vector{Symbol}, @view(expr.args[2:end]))
    # vecargs = [isa(arg, Symbol) ? :($arg::Union{T,Vec{W,T}})    : arg for arg ∈ args]
    mixargs = [isa(arg, Symbol) ? :($arg::Union{T,AbstractSIMDVector{W,T}})    : arg for arg ∈ args]
    svecargs = [isa(arg, Symbol) ? :($arg::SVec{W,T})    : arg for arg ∈ args]
    pt = [Symbol(:T_,i) for i ∈ eachindex(args)]
    promoteargs = [isa(args[i], Symbol) ? :($(args[i])::Union{$(pt[i]),AbstractSIMDVector{W,$(pt[i])}}) : args[i] for i ∈ eachindex(args)]
    q = quote
        @inline $f($(svecargs...)) where {W,T} = $(Expr(:call, :SVec, Expr(:call, f, [ Expr(:call, :extract_data, arg) for arg ∈ args  ]... ) ))
        # @inline $f($(mixargs...)) where {W,T} = $(Expr(:call, f, [ Expr(:call, :vbroadcast, Expr(:curly,:SVec,:W,:T), arg) for arg ∈ args  ]... ))
        @inline function $f($(promoteargs...)) where {W,$(pt...)}
            T = $(Expr(:call, :promote_type, pt...))
            $(Expr(:call, f, [ Expr(:call, :vconvert, Expr(:curly,:SVec,:W,:T), arg) for arg ∈ args  ]... ))
        end
    end
    esc(q)
#     esc(Expr(
#         :block,
#         Expr(:macrocall,# create all SVec definition
#              Symbol("@inline"), LineNumberNode(@__LINE__, @__FILE__),
#              Expr(:(=),
#                   Expr(:where, Expr(:call, f, svecargs...), :W, :T),
#                   Expr(:block,Expr(:call, :SVec, Expr(:call, f, [ Expr(:call, :extract_data, arg) for arg ∈ args  ]... ) ))
#                   )
#              ),
#         # Expr(:macrocall,# call all-Vec definition, promoting scalars
#              # Symbol("@inline"), LineNumberNode(@__LINE__, @__FILE__),
#              # Expr(:(=),
#                   # Expr(:where, Expr(:call, f, vecargs...), :W, :T),
#                   # Expr(:block,Expr(:call, f, [ Expr(:call, :vbroadcast, Expr(:curly,:Vec,:W,:T), arg) for arg ∈ args  ]... ))
#                   # )
#              # ),
#         Expr(:macrocall,# call all-SVec definition, promoting scalars and Vecs
#              Symbol("@inline"), LineNumberNode(@__LINE__, @__FILE__),
#              Expr(:(=),
#                   Expr(:where, Expr(:call, f, mixargs...), :W, :T),
#                   Expr(:block,Expr(:call, f, [ Expr(:call, :vbroadcast, Expr(:curly,:SVec,:W,:T), arg) for arg ∈ args  ]... ))
#                   )
#              )
#     ))
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
# include("arithmeticwithconsts.jl")
include("precompile.jl")
_precompile_()

function __init__()
    _precompile_()
end

end # module
