function _precompile_()
    ccall(:jl_generating_output, Cint, ()) == 1 || return nothing
    precompile(Tuple{typeof(SIMDPirates.capture_muladd),Expr,Nothing,Expr})
    precompile(Tuple{typeof(SIMDPirates.capture_muladd),Expr,Nothing,Symbol})
    precompile(Tuple{typeof(SIMDPirates.getneutral),Symbol,Type{Float32}})
    precompile(Tuple{typeof(SIMDPirates.getneutral),Symbol,Type{Int32}})
    precompile(Tuple{typeof(SIMDPirates.getneutral),Symbol,Type{Int64}})
    precompile(Tuple{typeof(SIMDPirates.llvmconst),Int,Type{T} where T,Int})
    precompile(Tuple{typeof(SIMDPirates.llvmins),Symbol,Int,Type{T} where T})
    precompile(Tuple{typeof(SIMDPirates.mulexpr),SubArray{Any,1,Array{Any,1},Tuple{UnitRange{Int}},true}})
    precompile(Tuple{typeof(SIMDPirates.subvector),String,Int,String,String,Int,Int})
    precompile(Tuple{typeof(SIMDPirates.vsub),VectorizationBase.SVec{8,Float32},VectorizationBase.SVec{8,Float64}})
end
