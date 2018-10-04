module SIMDPirates



export  Vec,
        broadcast,
        valloc,
        vload,
        vloada,
        vstore,
        vstorea,
        shufflevector,
        vifelse,
        vfma

include("type_definitions.jl")
include("llvm_utils.jl")
include("llvmwrap.jl")
include("conditionals.jl")
include("integer_arithmetic.jl")
include("floating_point_arithmetic.jl")
include("memory.jl")
include("shufflevector.jl")



end # module
