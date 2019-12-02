function _precompile_()
    ccall(:jl_generating_output, Cint, ()) == 1 || return nothing
    precompile(Tuple{typeof(SIMDPirates.llvmconst),Int64,Type,Float64})
    precompile(Tuple{typeof(SIMDPirates.llvmconst),Int64,Type,Int64})
end
