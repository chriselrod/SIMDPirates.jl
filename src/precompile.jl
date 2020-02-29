function _precompile_()
    ccall(:jl_generating_output, Cint, ()) == 1 || return nothing
    precompile(Tuple{typeof(SIMDPirates.getneutral),Symbol,Type{Float32}})
    precompile(Tuple{typeof(SIMDPirates.getneutral),Symbol,Type{Float64}})
    precompile(Tuple{typeof(SIMDPirates.getneutral),Symbol,Type{Int32}})
    precompile(Tuple{typeof(SIMDPirates.getneutral),Symbol,Type{Int64}})
    precompile(Tuple{typeof(SIMDPirates.llvmconst),Int64,Type{T} where T,Int64})
    precompile(Tuple{typeof(SIMDPirates.llvmins),Symbol,Int64,Type{T} where T})
    precompile(Tuple{typeof(SIMDPirates.mulexpr),SubArray{Any,1,Array{Any,1},Tuple{UnitRange{Int64}},true}})
    precompile(Tuple{typeof(SIMDPirates.subvector),String,Int64,String,String,Int64,Int64})
    precompile(Tuple{typeof(SIMDPirates.vsub),VectorizationBase.SVec{8,Float32},VectorizationBase.SVec{8,Float64}})
end
