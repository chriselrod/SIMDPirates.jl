module SIMDPirates

using MacroTools: postwalk, @capture

export  Vec, VE,
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

include("type_definitions.jl")
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
