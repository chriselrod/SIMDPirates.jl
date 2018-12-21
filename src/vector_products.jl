abstract type AbstractVectorProduct{N,T} <: AbstractStructVec{N,T} end
struct VecProduct{N,T} <: AbstractVectorProduct{N,T}
    v1::Vec{N,T}
    v2::Vec{N,T}
end
const VecOrProd{N,T} = Union{Vec{N,T},VecProduct{N,T}}
struct SVecProduct{N,T} <: AbstractVectorProduct{N,T}
    v1::Vec{N,T}
    v2::Vec{N,T}
end
@inline VectorizationBase.extract_data(v::AbstractVectorProduct) = llvmwrap(Val{:(*)}, v.v1, v.v2)
@inline SVec(x::VecProduct) = SVecProduct(x.v1, x.v2)
@inline SVec(x::SVecProduct) = x
