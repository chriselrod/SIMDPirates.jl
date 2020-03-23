using SIMDPirates
using Test



pkgdir(pkg::String) = abspath(joinpath(dirname(Base.find_package(pkg)), ".."))
include(joinpath(pkgdir("AccurateArithmetic"), "test", "runtests.jl"))


