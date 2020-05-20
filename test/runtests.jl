using SIMDPirates
using Test

pkgdir(pkg::String) = abspath(joinpath(dirname(Base.find_package(pkg)), ".."))

@testset "SIMDPirates.jl" begin
@test isempty(detect_unbound_args(SIMDPirates))
#pkgs = ["AccurateArithmetic", "MCMCChainSummaries", "VectorizedRNG"]
pkgs = ["AccurateArithmetic", "VectorizedRNG"]
for pkg ∈ pkgs
    include(joinpath(pkgdir(pkg), "test", "runtests.jl"))
end

sxi = SVec(Core.VecElement.((-1,0,1,2)))
@test sign(sxi) === SVec(Core.VecElement.((-1,0,1,1)))
@test sum(sxi) == 2
@test prod(sxi) == 0
@test maximum(sxi) == 2
@test minimum(sxi) == -1
@test SIMDPirates.addscalar(sxi, 2) === SVec(Core.VecElement.((1,0,1,2)))
@test SIMDPirates.mulscalar(sxi, 2) === SVec(Core.VecElement.((-2,0,1,2)))
@test SIMDPirates.maxscalar(sxi, 2) === SVec(Core.VecElement.((2,0,1,2)))
@test SIMDPirates.maxscalar(sxi, -2) === SVec(Core.VecElement.((-1,0,1,2)))
@test SIMDPirates.minscalar(sxi, 2) === SVec(Core.VecElement.((-1,0,1,2)))
@test SIMDPirates.minscalar(sxi, -2) === SVec(Core.VecElement.((-2,0,1,2)))

end

