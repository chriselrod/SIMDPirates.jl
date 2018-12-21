module SIMDPirates


using MacroTools: postwalk, @capture
using VectorizationBase

export  Vec, SVec, VE,
        @pirate,
        vbroadcast,
        valloc,
        vload,
        vloada,
        vstore,
        vstorea,
        shufflevector,
        vifelse,
        vfma,
        vmuladd


function vector_args(args)
    vecargs =       [isa(arg, Symbol) ? :($arg::VecOrProd{N,T})             : arg for arg ∈ args]
    abstractargs =  [isa(arg, Symbol) ? :($arg::AbstractSIMDVector{N,T})    : arg for arg ∈ args]
    structargs =    [isa(arg, Symbol) ? :($arg::AbstractStructVec{N,T})     : arg for arg ∈ args]
    vecargs, abstractargs, structargs
end

macro vectordef(rename, expr)
    if @capture(expr, (f_(args__) where {R__} = body_) | (function f_(args__) where {R__}; body_ end))
        vecargs, abstractargs, structargs = vector_args(args)
        q = quote
            @inline $rename($(vecargs...)) where {$(R...)} = $body
            @inline $rename($(abstractargs...)) where {$(R...)} = SVec($body)
            @inline $f($(structargs...)) where {$(R...)} = SVec($body)
        end
    else
        @show expr
        throw("Failed to match expression $expr")
    end
    esc(q)
end

include("type_definitions.jl")
include("vector_products.jl")
include("llvm_utils.jl")
include("llvmwrap.jl")
include("conditionals.jl")
include("integer_arithmetic.jl")
include("floating_point_arithmetic.jl")
include("memory.jl")
include("shufflevector.jl")
include("pirate.jl")
# include("restrict_simd.jl")

end # module
