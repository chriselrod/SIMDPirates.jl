function _precompile_()
    ccall(:jl_generating_output, Cint, ()) == 1 || return nothing
    precompile(Tuple{typeof(SIMDPirates.getneutral),Symbol,Type{Bool}})
    precompile(Tuple{typeof(SIMDPirates.getneutral),Symbol,Type{Float64}})
    precompile(Tuple{typeof(SIMDPirates.llvmconst),Int64,Type{T} where T,Float64})
    precompile(Tuple{typeof(SIMDPirates.llvmconst),Int64,Type{T} where T,Int64})
    precompile(Tuple{typeof(SIMDPirates.llvmins),Symbol,Int64,Type{T} where T})
    precompile(Tuple{typeof(SIMDPirates.scalar2vector),String,Int64,String,String})
    precompile(Tuple{typeof(SIMDPirates.subvector),String,Int64,String,String,Int64,Int64})
end
