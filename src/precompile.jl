function _precompile_()
    ccall(:jl_generating_output, Cint, ()) == 1 || return nothing
    precompile(Tuple{typeof(SIMDPirates.llvmconst),Int64,Type,Float64})
    precompile(Tuple{typeof(SIMDPirates.llvmconst),Int64,Type,Int64})
    precompile(Tuple{typeof(SIMDPirates.scalar2vector),String,Int64,String,String})
    precompile(Tuple{typeof(SIMDPirates.subvector),String,Int64,String,String,Int64,Int64})
end
