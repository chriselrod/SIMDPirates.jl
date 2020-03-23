using SIMDPirates
using Test

pkgdir(pkg::String) = abspath(joinpath(dirname(Base.find_package(pkg)), ".."))

pkgs = ["AccurateArithmetic", "MCMCChainSummaries", "VectorizedRNG"]
for pkg ∈ pkgs
    include(joinpath(pkgdir(pkg), "test", "runtests.jl"))
end


