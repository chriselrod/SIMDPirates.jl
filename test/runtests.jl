using SIMDPirates, VectorizationBase
using Test


for T in (Float16, Float32, Float64)
    @test sizeof(int_type(T)) == sizeof(T)
    @test sizeof(uint_type(T)) == sizeof(T)
    @test significand_bits(T) + exponent_bits(T) + sign_bits(T) == 8*sizeof(T)
    @test significand_mask(T) | exponent_mask(T) | sign_mask(T) ==
        typemax(uint_type(T))
    @test significand_mask(T) ⊻ exponent_mask(T) ⊻ sign_mask(T) ==
        typemax(uint_type(T))

    W = VectorizationBase.pick_vector_width(T)
    x = ntuple(W) do _ Core.VecElement(randn(T)) end
    y = ntuple(W) do _ Core.VecElement(randn(T)) end
    @test SIMDPirates.vall(
        SIMDPirates.vgreater_or_equal(x,y) === ntuple(_ -> Core.VecElement(x[w].value >= y[w].value), W)
    )
    @test SIMDPirates.vecbool_to_mask( SIMDPirates.vless(x,y) ) === SIMDPirates.vgreater_mask(y, x)
end
