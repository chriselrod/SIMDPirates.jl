using SIMDPirates
using Test

pkgdir(pkg::String) = abspath(joinpath(dirname(Base.find_package(pkg)), ".."))

#pkgs = ["AccurateArithmetic", "MCMCChainSummaries", "VectorizedRNG"]
pkgs = ["AccurateArithmetic", "VectorizedRNG"]
for pkg ∈ pkgs
    include(joinpath(pkgdir(pkg), "test", "runtests.jl"))
end


sx = SVec(Core.VecElement.((-1,0,1,2)))
@test sign(sx) === SVec(Core.VecElement.((-1,0,1,1)))
@test sum(sx) == 2
@test prod(sx) == 0
@test maximum(sx) == 2
@test minimum(sx) == -1


