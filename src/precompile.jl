function _precompile_()
    ccall(:jl_generating_output, Cint, ()) == 1 || return nothing
    precompile(Tuple{typeof(SIMDPirates.intrangetuple),Int64,Type{Float32}})
    precompile(Tuple{typeof(SIMDPirates.llvmconst),Int64,Type{Bool},Bool})
    precompile(Tuple{typeof(SIMDPirates.llvmconst),Int64,Type{T} where T,Int64})
    precompile(Tuple{typeof(SIMDPirates.mulexpr),SubArray{Any,1,Array{Any,1},Tuple{UnitRange{Int64}},true}})
    precompile(Tuple{typeof(SIMDPirates.sparse_strided_ptr_index),Core.SimpleVector,Int64,Type{Int64}})
end
