using SIMDPirates
using Base.Test


for T in (Float16, Float32, Float64)
    @test sizeof(int_type(T)) == sizeof(T)
    @test sizeof(uint_type(T)) == sizeof(T)
    @test significand_bits(T) + exponent_bits(T) + sign_bits(T) == 8*sizeof(T)
    @test significand_mask(T) | exponent_mask(T) | sign_mask(T) ==
        typemax(uint_type(T))
    @test significand_mask(T) ⊻ exponent_mask(T) ⊻ sign_mask(T) ==
        typemax(uint_type(T))
end
