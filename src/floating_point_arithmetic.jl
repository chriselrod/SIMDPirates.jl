
# Floating point arithmetic functions

for op ∈ (
        :(+), :(-),
        :abs, :floor, :ceil, :round,
        :sin, :cos,
        :exp, :exp2, :inv, :log, :log10, :log2,
        :sqrt, :trunc
    )
    rename = VECTOR_SYMBOLS[op]
    @eval begin
        @inline $rename(s1::FloatingTypes) = $op(s1)

        @vectordef $rename function Base.$op(v1) where {N,T<:FloatingTypes}
            llvmwrap(Val{$(QuoteNode(op))}, extract_data(v1))
        end

        # @inline $rename(v1::Vec{N,T}) where {N,T<:FloatingTypes} =
        #     llvmwrap(Val{$(QuoteNode(op))}, v1)
        # @inline $rename(v1::AbstractStructVec{N,T}) where {N,T<:FloatingTypes} =
        #     SVec(llvmwrap(Val{$(QuoteNode(op))}, extract_data(v1)))
        # @inline Base.$op(v1::AbstractStructVec{N,T}) where {N,T<:FloatingTypes} =
        #     SVec(llvmwrap(Val{$(QuoteNode(op))}, extract_data(v1)))
    end
end
@inline vexp10(s1::FloatingTypes) = exp10(s1)
@vectordef vexp10 function Base.exp10(v1) where {N,T <: FloatingTypes}
    vpow(vbroadcast(Vec{N,T}, 10), extract_data(v1))
end
# @inline vexp10(v1::AbstractSIMDVector{N,T}) where {N,T<:FloatingTypes} = vpow(vbroadcast(Vec{N,T}, 10), v1)
@inline vsign(s1::FloatingTypes) = sign(s1)
@vectordef vsign function Base.sign(v1) where {N,T<:FloatingTypes}
    vifelse(
        extract_data(v1) == vbroadcast(Vec{N,T}, zero(T)),
        vbroadcast(Vec{N,T}, zero(T)),
        copysign(vbroadcast(Vec{N,T}, one(T)), extract_data(v1))
    )
end
# @inline function vsign(v1::AbstractStructVec{N,T}) where {N,T<:FloatingTypes}
#     SVec(vsign(extract_data(v1)))
# end
# @inline function Base.sign(v1::AbstractStructVec{N,T}) where {N,T<:FloatingTypes}
#     SVec(vsign(extract_data(v1)))
# end

for op ∈ (
        :(-), :(/), :(%), :(^),
        :copysign, :max, :min
    )
    rename = VECTOR_SYMBOLS[op]
    @eval begin
        @inline $rename(s1::FloatingTypes, s2::FloatingTypes) = $op(s1, s2)

        @vectordef $rename function Base.$op(v1, v2) where {N,T <: FloatingTypes}
            llvmwrap(Val{$(QuoteNode(op))}, extract_data(v1), extract_data(v2))
        end

        # @inline $rename(v1::Vec{N,T}, v2::Vec{N,T}) where {N,T<:FloatingTypes} =
        #     llvmwrap(Val{$(QuoteNode(op))}, v1, v2)
        # @inline $rename(v1::AbstractSIMDVector{N,T}, v2::AbstractSIMDVector{N,T}) where {N,T<:FloatingTypes} =
        #     SVec(llvmwrap(Val{$(QuoteNode(op))}, extract_data(v1), extract_data(v2)))
        # @inline Base.$op(v1::AbstractSIMDVector{N,T}, v2::AbstractSIMDVector{N,T}) where {N,T<:FloatingTypes} =
        #     SVec(llvmwrap(Val{$(QuoteNode(op))}, extract_data(v1), extract_data(v2)))
    end
end

struct VecProduct{N,T} <: AbstractStructVec{N,T}
    v1::Vec{N,T}
    v2::Vec{N,T}
end
@inline extract_data(v::VecProduct) = llvmwrap(Val{:(*)}, v.v1, v.v2)

let op = :(*)
    rename = VECTOR_SYMBOLS[op]
    @eval begin
        @inline $rename(s1::FloatingTypes, s2::FloatingTypes) = $op(s1, s2)
        # @inline $rename(v1::Vec{N,T}, v2::Vec{N,T}) where {N,T<:FloatingTypes} =
        #     VecProduct(v1, v2)
        @inline $rename(v1::AbstractSIMDVector{N,T}, v2::AbstractSIMDVector{N,T}) where {N,T<:FloatingTypes} =
            VecProduct(extract_data(v1), extract_data(v2))
        @inline Base.$op(v1::AbstractStructVec{N,T}, v2::AbstractStructVec{N,T}) where {N,T<:FloatingTypes} =
            VecProduct(extract_data(v1), extract_data(v2))
    end
end

let op = :(+)
    rename = VECTOR_SYMBOLS[op]
    @eval begin
        @inline $rename(s1::FloatingTypes, s2::FloatingTypes) = $op(s1, s2)
        @inline $rename(v1::Vec{N,T}, v2::Vec{N,T}) where {N,T<:FloatingTypes} =
            llvmwrap(Val{$(QuoteNode(op))}, v1, v2)
        @inline $rename(v1::AbstractSIMDVector{N,T}, v2::AbstractSIMDVector{N,T}) where {N,T<:FloatingTypes} =
            SVec(llvmwrap(Val{$(QuoteNode(op))}, extract_data(v1), extract_data(v2)))
        @inline Base.$op(v1::AbstractSIMDVector{N,T}, v2::AbstractSIMDVector{N,T}) where {N,T<:FloatingTypes} =
            SVec(llvmwrap(Val{$(QuoteNode(op))}, extract_data(v1), extract_data(v2)))

        @inline $rename(v::VecProduct{N,T}, v3::Vec{N,T}) where {N,T<:FloatingTypes} =
            llvmwrap(Val{:fma}, v.v1, v.v2, v3)
        @inline $rename(v3::Vec{N,T}, v::VecProduct{N,T}) where {N,T<:FloatingTypes} =
            llvmwrap(Val{:fma}, v.v1, v.v2, v3)
        @inline $rename(v::VecProduct{N,T}, v3::AbstractStructVec{N,T}) where {N,T<:FloatingTypes} =
            SVec(llvmwrap(Val{:fma}, v.v1, v.v2, extract_data(v3)))
        @inline $rename(v3::AbstractStructVec{N,T}, v::VecProduct{N,T}) where {N,T<:FloatingTypes} =
            SVec(llvmwrap(Val{:fma}, v.v1, v.v2, extract_data(v3)))

        @inline Base.$op(v::VecProduct{N,T}, v3::Vec{N,T}) where {N,T<:FloatingTypes} =
            llvmwrap(Val{:fma}, v.v1, v.v2, v3)
        @inline Base.$op(v3::Vec{N,T}, v::VecProduct{N,T}) where {N,T<:FloatingTypes} =
            llvmwrap(Val{:fma}, v.v1, v.v2, v3)
        @inline Base.$op(v::VecProduct{N,T}, v3::AbstractStructVec{N,T}) where {N,T<:FloatingTypes} =
            SVec(llvmwrap(Val{:fma}, v.v1, v.v2, extract_data(v3)))
        @inline Base.$op(v3::AbstractStructVec{N,T}, v::VecProduct{N,T}) where {N,T<:FloatingTypes} =
            SVec(llvmwrap(Val{:fma}, v.v1, v.v2, extract_data(v3)))


        @inline $rename(vp1::VecProduct{N,T}, vp2::VecProduct{N,T}) where {N,T<:FloatingTypes} =
            SVec(llvmwrap(Val{:fma}, vp1.v1, vp1.v2, extract_data(vp2)))
        @inline Base.$op(vp1::VecProduct{N,T}, vp2::VecProduct{N,T}) where {N,T<:FloatingTypes} =
            SVec(llvmwrap(Val{:fma}, vp1.v1, vp1.v2, extract_data(vp2)))
    end
end

@inline vpow(s1::FloatingTypes, x2::Integer) = s1^x2
@vectordef vpow function Base.:^(v1, x2::Integer) where {N,T<:FloatingTypes}
    llvmwrap(Val{:powi}, extract_data(v1), Int(x2))
end
# @inline function vpow(v1::Vec{N,T}, x2::Integer) where {N,T<:FloatingTypes}
#     llvmwrap(Val{:powi}, v1, Int(x2))
# end
# @inline function vpow(v1::AbstractStructVec{N,T}, x2::Integer) where {N,T<:FloatingTypes}
#     SVec(llvmwrap(Val{:powi}, extract_data(v1), Int(x2)))
# end
# @inline function Base.:^(v1::AbstractStructVec{N,T}, x2::Integer) where {N,T<:FloatingTypes}
#     SVec(llvmwrap(Val{:powi}, extract_data(v1), Int(x2)))
# end

@inline vflipsign(s1::FloatingTypes, s2::FloatingTypes) = flipsign(s1, s2)

@vectordef vflipsign function Base.flipsign(v1, v2) where {N,T<:FloatingTypes}
    vifelse(vsignbit(extract_data(v2)), vsub(extract_data(v1)), extract_data(v1))
end

# @inline vflipsign(v1::Vec{N,T}, v2::Vec{N,T}) where {N,T<:FloatingTypes} =
#     vifelse(vsignbit(v2), -v1, v1)
# @inline vflipsign(v1::AbstractSIMDVector{N,T}, v2::AbstractSIMDVector{N,T}) where {N,T<:FloatingTypes} =
#     SVec(vifelse(vsignbit(v2), -v1, v1))
# @inline Base.flipsign(v1::AbstractStructVec{N,T}, v2::AbstractStructVec{N,T}) where {N,T<:FloatingTypes} =
#     SVec(vifelse(vsignbit(v2), -v1, v1))

for op ∈ (:fma, :muladd)
    rename = VECTOR_SYMBOLS[op]
    @eval begin

        @vectordef $rename function Base.$op(v1, v2, v3) where {N,T<:FloatingTypes}
            llvmwrap(Val{$(QuoteNode(op))}, extract_data(v1), extract_data(v2), extract_data(v3))
        end

        # scalar default already set in integer_arithmetic.jl
        # @inline function $rename(v1::Vec{N,T},
        #         v2::Vec{N,T}, v3::Vec{N,T}) where {N,T<:FloatingTypes}
        #     llvmwrap(Val{$(QuoteNode(op))}, v1, v2, v3)
        # end
        # @inline function $rename(v1::AbstractSIMDVector{N,T},
        #         v2::AbstractSIMDVector{N,T}, v3::AbstractSIMDVector{N,T}) where {N,T<:FloatingTypes}
        #     SVec(llvmwrap(Val{$(QuoteNode(op))}, extract_data(v1), extract_data(v2), extract_data(v3)))
        # end
        # @inline function Base.$op(v1::AbstractStructVec{N,T},
        #         v2::AbstractStructVec{N,T}, v3::AbstractStructVec{N,T}) where {N,T<:FloatingTypes}
        #     SVec(llvmwrap(Val{$(QuoteNode(op))}, extract_data(v1), extract_data(v2), extract_data(v3)))
        # end
    end
end

# Type promotion



# Promote scalars of all ScalarTypes to vectors of FloatingTypes, leaving the
# vector type unchanged

for op ∈ (
        :(==), :(!=), :(<), :(<=), :(>), :(>=),
        :+, :-, :*, :/, :^,
        :copysign, :flipsign, :max, :min, :%
    )
    rename = VECTOR_SYMBOLS[op]
    @eval begin
        @vectordef $rename function Base.$op(s1::ScalarTypes, v2) where {N,T<:FloatingTypes}
            $rename(vbroadcast(Vec{N,T}, s1), extract_data(v2))
        end
        @vectordef $rename function Base.$op(v1, s2::ScalarTypes) where {N,T<:FloatingTypes}
            $rename(extract_data(v1), vbroadcast(Vec{N,T}, s2))
        end

        # @inline $rename(s1::ScalarTypes, v2::Vec{N,T}) where {N,T<:FloatingTypes} =
        #     $rename(vbroadcast(Vec{N,T}, s1), v2)
        # @inline $rename(v1::Vec{N,T}, s2::ScalarTypes) where {N,T<:FloatingTypes} =
        #     $rename(v1, vbroadcast(Vec{N,T}, s2))
    end
end
# @vectordef vifelse function Base.ifelse(c::AbstractSIMDVector{N,Bool}, s1::ScalarTypes, v2) where {N,T<:FloatingTypes}
#     vifelse(extract_data(c), vbroadcast(Vec{N,T}, s1), extract_data(v2))
# end
# @vectordef vifelse function Base.ifelse(c::AbstractSIMDVector{N,Bool}, v1, s2::ScalarTypes) where {N,T<:FloatingTypes}
#     vifelse(extract_data(c), extract_data(v1), vbroadcast(Vec{N,T}, s2))
# end
@inline function vifelse(c::Vec{N,Bool}, s1::ScalarTypes, v2::Vec{N,T}) where {N,T<:FloatingTypes}
    vifelse(c, vbroadcast(Vec{N,T}, s1), v2)
end
@inline function vifelse(c::Vec{N,Bool}, v1::Vec{N,T}, s2::ScalarTypes) where {N,T<:FloatingTypes}
    vifelse(c, v1, vbroadcast(Vec{N,T}, s2))
end
@inline function vifelse(c::AbstractSIMDVector{N,Bool}, s1::ScalarTypes, v2::AbstractSIMDVector{N,T}) where {N,T<:FloatingTypes}
    SVec(vifelse(extract_data(c), vbroadcast(Vec{N,T}, s1), extract_data(v2)))
end
@inline function vifelse(c::AbstractSIMDVector{N,Bool}, v1::AbstractSIMDVector{N,T}, s2::ScalarTypes) where {N,T<:FloatingTypes}
    SVec(vifelse(extract_data(c), extract_data(v1), vbroadcast(Vec{N,T}, s2)))
end
# @inline vifelse(c::Vec{N,Bool}, s1::ScalarTypes,
#         v2::Vec{N,T}) where {N,T<:FloatingTypes} =
#     vifelse(c, vbroadcast(Vec{N,T}, s1), v2)
# @inline vifelse(c::Vec{N,Bool}, v1::Vec{N,T},
#         s2::ScalarTypes) where {N,T<:FloatingTypes} =
# vifelse(c, v1, vbroadcast(Vec{N,T}, s2))

for op ∈ (:fma, :muladd)
    rename = VECTOR_SYMBOLS[op]
    @eval begin
        @vectordef $rename function Base.$op(s1::ScalarTypes, v2, v3) where {N,T}
            $rename(vbroadcast(Vec{N,T}, s1), extract_data(v2), extract_data(v3))
        end
        @vectordef $rename function Base.$op(v1, s2::ScalarTypes, v3) where {N,T}
            $rename(extract_data(v1), vbroadcast(Vec{N,T}, s2), extract_data(v3))
        end
        @vectordef $rename function Base.$op(v1, v2, s3::ScalarTypes) where {N,T}
            $rename(extract_data(v1), extract_data(v2), vbroadcast(Vec{N,T}, s3))
        end


        @vectordef $rename function Base.$op(s1::ScalarTypes, s2::ScalarTypes, v3) where {N,T}
            $rename(vbroadcast(Vec{N,T}, s1), vbroadcast(Vec{N,T}, s2), extract_data(v3))
        end
        @vectordef $rename function Base.$op(s1::ScalarTypes, v2, s3::ScalarTypes) where {N,T}
            $rename(vbroadcast(Vec{N,T}, s1), extract_data(v2), vbroadcast(Vec{N,T}, s3))
        end
        @vectordef $rename function Base.$op(v1, s2::ScalarTypes, s3::ScalarTypes) where {N,T}
            $rename(extract_data(v1), vbroadcast(Vec{N,T}, s2), vbroadcast(Vec{N,T}, s3))
        end

        # @inline $rename(s1::ScalarTypes, v2::Vec{N,T},
        #         v3::Vec{N,T}) where {N,T<:FloatingTypes} =
        #     $rename(vbroadcast(Vec{N,T}, s1), v2, v3)
        # @inline $rename(v1::Vec{N,T}, s2::ScalarTypes,
        #         v3::Vec{N,T}) where {N,T<:FloatingTypes} =
        #     $rename(v1, vbroadcast(Vec{N,T}, s2), v3)
        # @inline $rename(s1::ScalarTypes, s2::ScalarTypes,
        #         v3::Vec{N,T}) where {N,T<:FloatingTypes} =
        #     $rename(vbroadcast(Vec{N,T}, s1), vbroadcast(Vec{N,T}, s2), v3)
        # @inline $rename(v1::Vec{N,T}, v2::Vec{N,T},
        #         s3::ScalarTypes) where {N,T<:FloatingTypes} =
        #     $rename(v1, v2, vbroadcast(Vec{N,T}, s3))
        # @inline $rename(s1::ScalarTypes, v2::Vec{N,T},
        #         s3::ScalarTypes) where {N,T<:FloatingTypes} =
        #     $rename(vbroadcast(Vec{N,T}, s1), v2, vbroadcast(Vec{N,T}, s3))
        # @inline $rename(v1::Vec{N,T}, s2::ScalarTypes,
        #         s3::ScalarTypes) where {N,T<:FloatingTypes} =
        #     $rename(v1, vbroadcast(Vec{N,T}, s2), vbroadcast(Vec{N,T}, s3))
    end
end

# Reduction operations

# TODO: map, mapreduce

function getneutral(op::Symbol, ::Type{T}) where T
    zs = Dict{Symbol,T}()
    if T <: IntegerTypes
        zs[:&] = ~T(0)
        zs[:|] = T(0)
    end
    zs[:max] = typemin(T)
    zs[:min] = typemax(T)
    zs[:+] = T(0)
    zs[:*] = T(1)
    zs[op]
end

nextpow2(n) = nextpow(2, n)

# We cannot pass in the neutral element via Val{}; if we try, Julia refuses to
# inline this function, which is then disastrous for performance
@generated function llvmwrapreduce(::Type{Val{Op}}, v::Vec{N,T}) where {Op,N,T}
    @assert isa(Op, Symbol)
    z = getneutral(Op, T)
    typ = llvmtype(T)
    decls = []
    instrs = []
    n = N
    nam = "%0"
    nold,n = n,nextpow2(n)
    if n > nold
        namold,nam = nam,"%vec_$n"
        append!(instrs,
            extendvector(namold, nold, typ, n, n-nold,
                llvmtypedconst(T,z), nam))
    end
    while n > 1
        nold,n = n, div(n, 2)
        namold,nam = nam,"%vec_$n"
        vtyp = "<$n x $typ>"
        ins = llvmins(Val{Op}, n, T)
        append!(instrs, subvector(namold, nold, typ, "$(nam)_1", n, 0))
        append!(instrs, subvector(namold, nold, typ, "$(nam)_2", n, n))
        if ins[1] == '@'
            push!(decls, "declare $vtyp $ins($vtyp, $vtyp)")
            push!(instrs,
                "$nam = call $vtyp $ins($vtyp $(nam)_1, $vtyp $(nam)_2)")
        else
            push!(instrs, "$nam = $ins $vtyp $(nam)_1, $(nam)_2")
        end
    end
    push!(instrs, "%res = extractelement <$n x $typ> $nam, i32 0")
    push!(instrs, "ret $typ %res")
    quote
        $(Expr(:meta, :inline))
        Base.llvmcall($((join(decls, "\n"), join(instrs, "\n"))),
            T, Tuple{Vec{N,T}}, v)
    end
end

@vectordef vall Base.all(v) where {N,T<:IntegerTypes} = llvmwrapreduce(Val{:&}, extract_data(v))
@vectordef vany Base.any(v) where {N,T<:IntegerTypes} = llvmwrapreduce(Val{:|}, extract_data(v))
@vectordef vmaximum function Base.maximum(v) where {N,T<:FloatingTypes}
    llvmwrapreduce(Val{:max}, extract_data(v))
end
@vectordef vminimum function Base.minimum(v) where {N,T<:FloatingTypes}
    llvmwrapreduce(Val{:min}, extract_data(v))
end
@vectordef vsum function Base.sum(v) where {N,T}
    llvmwrapreduce(Val{:max}, extract_data(v))
end
@vectordef vprod function Base.prod(v) where {N,T}
    llvmwrapreduce(Val{:min}, extract_data(v))
end

# @inline vall(v::Vec{N,T}) where {N,T<:IntegerTypes} = llvmwrapreduce(Val{:&}, v)
# @inline vany(v::Vec{N,T}) where {N,T<:IntegerTypes} = llvmwrapreduce(Val{:|}, v)
# @inline vmaximum(v::Vec{N,T}) where {N,T<:FloatingTypes} =
#     llvmwrapreduce(Val{:max}, v)
# @inline vminimum(v::Vec{N,T}) where {N,T<:FloatingTypes} =
#     llvmwrapreduce(Val{:min}, v)
# @inline vprod(v::Vec{N,T}) where {N,T} = llvmwrapreduce(Val{:*}, v)
# @inline vsum(v::Vec{N,T}) where {N,T} = llvmwrapreduce(Val{:+}, v)

@generated function vreduce(::Type{Val{Op}}, v::Vec{N,T}) where {Op,N,T}
    @assert isa(Op, Symbol)
    z = getneutral(Op, T)
    stmts = []
    n = N
    push!(stmts, :($(Symbol(:v,n)) = v))
    nold,n = n,nextpow2(n)
    if n > nold
        push!(stmts,
            :($(Symbol(:v,n)) = $(Expr(:tuple,
                [:($(Symbol(:v,nold))[$i]) for i in 1:nold]...,
                [:(VE($z)) for i in nold+1:n]...))))
    end
    while n > 1
        nold,n = n, div(n, 2)
        push!(stmts,
            :($(Symbol(:v,n,"lo")) = $(Expr(:tuple,
                [:($(Symbol(:v,nold))[$i]) for i in 1:n]...,))))
        push!(stmts,
            :($(Symbol(:v,n,"hi")) = $(Expr(:tuple,
                [:($(Symbol(:v,nold))[$i]) for i in n+1:nold]...))))
        push!(stmts,
            :($(Symbol(:v,n)) =
                $Op($(Symbol(:v,n,"lo")), $(Symbol(:v,n,"hi")))))
    end
    push!(stmts, :(v1[1].value))
    Expr(:block, Expr(:meta, :inline), stmts...)
end

@vectordef vmaximum function Base.maximum(v) where {N,T<:IntegerTypes}
    vreduce(Val{:max}, extract_data(v))
end
@vectordef vminimum function Base.minimum(v) where {N,T<:IntegerTypes}
    vreduce(Val{:min}, extract_data(v))
end

# @inline vmaximum(v::Vec{N,T}) where {N,T<:IntegerTypes} = vreduce(Val{:max}, v)
# @inline vminimum(v::Vec{N,T}) where {N,T<:IntegerTypes} = vreduce(Val{:min}, v)


@inline vmul(x,y,z...) = vmul(x,vmul(y,z...))
@inline vadd(x,y,z...) = vadd(x,vadd(y,z...))
