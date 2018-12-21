
# Integer arithmetic functions

for op ∈ (:(~), :(+), :(-))
    rename = VECTOR_SYMBOLS[op]
    @eval begin
        @inline $rename(s1::IntegerTypes) = $op($IntegerTypes)
        @vectordef $rename function Base.$op(v1) where {N,T<:IntegerTypes}
            llvmwrap(Val{$(QuoteNode(op))}, extract_data(v1))
        end
    end
end
# @inline vnot(s1::Bool) = !s1
@inline vnot(v1::Vec{N,Bool}) where {N} = vbitwise_not(v1)
@inline vnot(v1::AbstractStructVec{N,Bool}) where {N} = SVec(vbitwise_not(extract_data(v1)))
@inline Base.:!(v1::AbstractStructVec{N,Bool}) where {N} = SVec(vbitwise_not(extract_data(v1)))

@inline vabs(s1::IntTypes) = abs(s1)
@inline function vabs(v1::Vec{N,T}) where {N,T<:IntTypes}
    # s = -vbroadcast(Vec{N,T}, signbit(v1))
    s = vright_bitshift(v1, Val{8*sizeof(T)}())
    # Note: -v1 == ~v1 + 1
    vsub(vxor(s, v1), s)
end
@inline vabs(v1::AbstractStructVec{N,T}) where {N,T<:IntTypes} = SVec(vabs(extract_data(v1)))
@inline Base.abs(v1::AbstractStructVec{N,T}) where {N,T<:IntTypes} = SVec(vabs(extract_data(v1)))

@inline vabs(v1::AbstractSIMDVector{N,T}) where {N,T<:UIntTypes} = v1
@inline Base.abs(v1::AbstractStructVec{N,T}) where {N,T<:UIntTypes} = v1

# TODO: Try T(v1>0) - T(v1<0)
#       use a shift for v1<0
#       evaluate v1>0 as -v1<0 ?

@inline vsign(s1::IntTypes) = sign(s1)
@inline vsign(v1::Vec{N,T}) where {N,T<:IntTypes} =
    vifelse(vequal(v1, vbroadcast(Vec{N,T},0)), vbroadcast(Vec{N,T},0),
        vifelse(vless(v1, vbroadcast(Vec{N,T},0)), vbroadcast(Vec{N,T},-1), vbroadcast(Vec{N,T},1)))
@inline vsign(v1::Vec{N,T}) where {N,T<:UIntTypes} =
    vifelse(vequal(v1, vbroadcast(Vec{N,T},0)), vbroadcast(Vec{N,T},0), vbroadcast(Vec{N,T},1))
    @inline vsign(v1::AbstractStructVec) = SVec(vsign(extract_data(v1)))
    @inline Base.sign(v1::AbstractStructVec) = SVec(vsign(extract_data(v1)))


@inline vsignbit(s1::IntegerTypes) = signbit(s1)
@inline vsignbit(v1::Vec{N,T}) where {N,T<:IntTypes} = vless(v1, vbroadcast(Vec{N,T}, 0))
@inline vsignbit(v1::AbstractStructVec{N,T}) where {N,T<:IntTypes} = SVec(vsignbit(extract_data(v1)))
@inline Base.signbit(v1::AbstractStructVec{N,T}) where {N,T<:IntTypes} = SVec(vsignbit(extract_data(v1)))


@inline vsignbit(v1::Vec{N,T}) where {N,T<:UIntTypes} = vbroadcast(Vec{N,Bool}, false)
@inline vsignbit(v1::AbstractStructVec{N,T}) where {N,T<:UIntTypes} = vbroadcast(SVec{N,Bool}, false)
@inline Base.signbit(v1::AbstractStructVec{N,T}) where {N,T<:UIntTypes} = vbroadcast(SVec{N,Bool}, false)

for op ∈ (:(&), :(|), :(⊻), :(+), :(-), :(*), :(÷), :(%) )
    rename = VECTOR_SYMBOLS[op]
    @eval begin
        @inline $rename(s1::IntegerTypes, s2::IntegerTypes) = $op(s1, s2)

        @vectordef $rename function Base.$op(v1, v2) where {N,T<:IntegerTypes}
            llvmwrap(Val{$(QuoteNode(op))}, extract_data(v1), extract_data(v2))
        end

        # @inline $rename(v1::Vec{N,T}, v2::Vec{N,T}) where {N,T<:IntegerTypes} =
        #     llvmwrap(Val{$(QuoteNode(op))}, v1, v2)
    end
end

let op = :(%)
    rename = VECTOR_SYMBOLS[op]

    @eval @vectordef $rename function Base.$op(v1, ::Type{I}) where {N, T <: IntegerTypes, I <: IntegerTypes}
        ntuple(Val(N)) do n
            VE(v1[n].value % I)
        end
    end

end


@inline vcopysign(s1::IntegerTypes, s2::IntegerTypes) = copysign(s1, s2)
@inline vcopysign(v1::Vec{N,T}, v2::Vec{N,T}) where {N,T<:IntTypes} =
    vifelse(vsignbit(v2), -abs(v1), abs(v1))
@inline function vcopysign(v1::AbstractStructVec{N,T}, v2::AbstractStructVec{N,T}) where {N,T}
    SVec(vcopysign(extract_data(v1), extract_data(v2)))
end
@inline function Base.copysign(v1::AbstractStructVec{N,T}, v2::AbstractStructVec{N,T}) where {N,T}
    SVec(vcopysign(extract_data(v1), extract_data(v2)))
end
@inline vcopysign(v1::AbstractSIMDVector{N,T}, v2::AbstractSIMDVector{N,T}) where {N,T<:UIntTypes} = v1
@inline Base.copysign(v1::AbstractStructVec{N,T}, v2::AbstractStructVec{N,T}) where {N,T<:UIntTypes} = v1

@inline vflipsign(s1::IntegerTypes, s2::IntegerTypes) = flipsign(s1, s2)
@inline vflipsign(v1::Vec{N,T}, v2::Vec{N,T}) where {N,T<:IntTypes} =
    vifelse(vsignbit(v2), -v1, v1)
@inline function vflipsign(v1::AbstractStructVec{N,T}, v2::AbstractStructVec{N,T}) where {N,T<:IntTypes}
    SVec(vflipsign(extract_data(v1), extract_data(V2)))
end
@inline function Base.flipsign(v1::AbstractStructVec{N,T}, v2::AbstractStructVec{N,T}) where {N,T<:IntTypes}
    SVec(vflipsign(extract_data(v1), extract_data(V2)))
end

@inline vflipsign(v1::AbstractSIMDVector{N,T}, v2::AbstractSIMDVector{N,T}) where {N,T<:UIntTypes} = v1
@inline Base.flipsign(v1::AbstractStructVec{N,T}, v2::AbstractStructVec{N,T}) where {N,T<:UIntTypes} = v1

@inline vmax(s1::IntegerTypes, s2::IntegerTypes) = vmax(s1, s2)
@inline vmax(v1::Vec{N,T}, v2::Vec{N,T}) where {N,T<:IntegerTypes} =
    vifelse(v1>=v2, v1, v2)
@inline vmin(v1::Vec{N,T}, v2::Vec{N,T}) where {N,T<:IntegerTypes} =
    vifelse(v1>=v2, v2, v1)

vmuladd(s1::ScalarTypes, s2::ScalarTypes, s3::ScalarTypes) = muladd(s1,s2,s3)
@inline function vmuladd(v1::Vec{N,T}, v2::Vec{N,T}, v3::Vec{N,T}) where {N,T<:IntegerTypes}
    vadd(vmul(v1,v2),v3)
end
@inline function vmuladd(v1::AbstractSIMDVector{N,T}, v2::AbstractSIMDVector{N,T},
        v3::AbstractSIMDVector{N,T}) where {N,T<:IntegerTypes}
    SVec(vadd(vmul(extract_data(v1),extract_data(v2)),extract_data(v3)))
end
@inline function Base.muladd(v1::AbstractStructVec{N,T}, v2::AbstractStructVec{N,T},
        v3::AbstractStructVec{N,T}) where {N,T<:IntegerTypes}
    SVec(vadd(vmul(extract_data(v1),extract_data(v2)),extract_data(v3)))
end

vfma(s1::ScalarTypes, s2::ScalarTypes, s3::ScalarTypes) = fma(s1,s2,s3)
@inline function vfma(v1::Vec{N,T}, v2::Vec{N,T}, v3::Vec{N,T}) where {N,T<:IntegerTypes}
    vadd(vmul(v1,v2),v3)
end
@inline function vfma(v1::AbstractSIMDVector{N,T}, v2::AbstractSIMDVector{N,T},
        v3::AbstractSIMDVector{N,T}) where {N,T<:IntegerTypes}
    SVec(vadd(vmul(extract_data(v1),extract_data(v2)),extract_data(v3)))
end
@inline function Base.fma(v1::AbstractStructVec{N,T}, v2::AbstractStructVec{N,T},
        v3::AbstractStructVec{N,T}) where {N,T<:IntegerTypes}
    SVec(vadd(vmul(extract_data(v1),extract_data(v2)),extract_data(v3)))
end

# TODO: Handle negative shift counts
#       use vifelse
#       ensure vifelse is efficient
for op ∈ (:(<<), :(>>), :(>>>))
    rename = VECTOR_SYMBOLS[op]
    @eval begin
        @inline $rename(s1::IntegerTypes, s2::IntegerTypes) = $op(s1, s2)

        @vectordef $rename function Base.$op(v1, ::Val{I}) where {N,T<:IntegerTypes,I}
            llvmwrapshift(Val{$(QuoteNode(op))}, extract_data(v1), Val{I})
        end
        @vectordef $rename function Base.$op(v1, x2::Unsigned) where {N,T<:IntegerTypes}
            llvmwrapshift(Val{$(QuoteNode(op))}, extract_data(v1), x2)
        end
        @vectordef $rename function Base.$op(v1, x2::Int) where {N,T<:IntegerTypes}
            llvmwrapshift(Val{$(QuoteNode(op))}, extract_data(v1), x2)
        end
        @vectordef $rename function Base.$op(v1, x2::Integer) where {N,T<:IntegerTypes}
            llvmwrapshift(Val{$(QuoteNode(op))}, extract_data(v1), x2)
        end
        @vectordef $rename function Base.$op(v1, x2::Vec{N,U}) where {N,T<:IntegerTypes,U<:UIntTypes}
            llvmwrapshift(Val{$(QuoteNode(op))}, extract_data(v1), extract_data(v2))
        end
        @vectordef $rename function Base.$op(x1::T, v2) where {N,T<:IntegerTypes}
            $rename(vbroadcast(Vec{N,T}, x1), extract_data(v2))
        end


        # @inline $rename(v1::Vec{N,T}, ::Val{I}) where {N,T<:IntegerTypes,I} =
            # llvmwrapshift(Val{$(QuoteNode(op))}, v1, Val{I})
        # @inline $rename(v1::Vec{N,T}, x2::Unsigned) where {N,T<:IntegerTypes} =
        #     llvmwrapshift(Val{$(QuoteNode(op))}, v1, x2)
        # @inline $rename(v1::Vec{N,T}, x2::Int) where {N,T<:IntegerTypes} =
        #     llvmwrapshift(Val{$(QuoteNode(op))}, v1, x2)
        # @inline $rename(v1::Vec{N,T}, x2::Integer) where {N,T<:IntegerTypes} =
        #     llvmwrapshift(Val{$(QuoteNode(op))}, v1, x2)
        # @inline $rename(v1::Vec{N,T},
        #                  v2::Vec{N,U}) where {N,T<:IntegerTypes,U<:UIntTypes} =
        #     llvmwrapshift(Val{$(QuoteNode(op))}, v1, v2)
        # @inline $rename(v1::Vec{N,T},
        #                  v2::Vec{N,U}) where {N,T<:IntegerTypes,U<:IntegerTypes} =
        #     llvmwrapshift(Val{$(QuoteNode(op))}, v1, v2)
        # @inline $rename(x1::T, v2::Vec{N,T}) where {N,T<:IntegerTypes} =
        #     $rename(vbroadcast(Vec{N,T}, x1), v2)
    end
end


# Promote scalars of all IntegerTypes to vectors of IntegerTypes, leaving the
# vector type unchanged

for op ∈ (
        :(==), :(!=), :(<), :(<=), :(>), :(>=),
        :(&), :(|), :(⊻), :(+), :(-), :(*),
        :copysign, :div, :flipsign, :max, :min, :rem
    )
    rename = VECTOR_SYMBOLS[op]
    @eval begin
        @inline $rename(s1::ScalarTypes, s2::ScalarTypes) = $op(s1, s2)

        @vectordef $rename function Base.$op(s1::T, v2) where {N, T <: Union{Bool,IntegerTypes}}
            $rename(vbroadcast(Vec{N,T}, s1), extract_data(v2))
        end
        @vectordef $rename function Base.$op(v2, s1::T) where {N, T <: Union{Bool,IntegerTypes}}
            $rename(vbroadcast(Vec{N,T}, s1), extract_data(v2))
        end


        # @vectordef $rename function Base.$op(s1::Bool, v2) where {N, T <: Bool}
        #     $rename(vbroadcast(Vec{N,Bool},s1), extract_data(v2))
        # end
        # @vectordef $rename function Base.$op(v2, s1::Bool) where {N, T <: Bool}
        #     $rename(vbroadcast(Vec{N,Bool},s1), extract_data(v2))
        # end
        # @inline $rename(s1::Bool, v2::Vec{N,Bool}) where {N} =
        #     $op(vbroadcast(Vec{N,Bool},s1), v2)
        # @inline $rename(s1::IntegerTypes, v2::Vec{N,T}) where {N,T<:IntegerTypes} =
        #     $op(vbroadcast(Vec{N,T}, s1), v2)
        # @inline $rename(v1::Vec{N,T}, s2::IntegerTypes) where {N,T<:IntegerTypes} =
        #     $op(v1, vbroadcast(Vec{N,T}, s2))
        #
        # @inline $rename(s1::Bool, v2::AbstractStructVec{N,Bool}) where {N} =
        #     SVec($op(Vec{N,Bool}(s1), extract_data(v2)))
        # @inline $rename(s1::IntegerTypes, v2::AbstractStructVec{N,T}) where {N,T<:IntegerTypes} =
        #     SVec($op(vbroadcast(Vec{N,T}, s1), extract_data(v2)))
        # @inline $rename(v1::AbstractStructVec{N,T}, s2::IntegerTypes) where {N,T<:IntegerTypes} =
        #     SVec($op(extract_data(v1), vbroadcast(Vec{N,T}, s2)))
    end
end
# @vectordef vifelse function Base.ifelse(c::Vec{N,Bool}, s1::T, v2) where {N,T<:IntegerTypes}
#     vifelse(c, vbroadcast(Vec{N,T}, s1), extract_data(v2))
# end
# @vectordef vifelse function Base.ifelse(c::Vec{N,Bool}, v1, s2::T) where {N,T<:IntegerTypes}
#     vifelse(c, extract_data(v1), vbroadcast(Vec{N,T}, s2))
# end
@inline function vifelse(c::Vec{N,Bool}, s1::T, v2::Vec{N,T}) where {N,T<:IntegerTypes}
    vifelse(c, vbroadcast(Vec{N,T}, s1), v2)
end
@inline function vifelse(c::Vec{N,Bool}, v1::Vec{N,T}, s2::T) where {N,T<:IntegerTypes}
    vifelse(c, v1, vbroadcast(Vec{N,T}, s2))
end
@inline function vifelse(c::AbstractSIMDVector{N,Bool}, s1::T, v2::AbstractSIMDVector{N,T}) where {N,T<:IntegerTypes}
    SVec(vifelse(extract_data(c), vbroadcast(Vec{N,T}, s1), extract_data(v2)))
end
@inline function vifelse(c::AbstractSIMDVector{N,Bool}, v1::AbstractSIMDVector{N,T}, s2::T) where {N,T<:IntegerTypes}
    SVec(vifelse(extract_data(c), extract_data(v1), vbroadcast(Vec{N,T}, s2)))
end


# @inline function vifelse(c::Vec{N,Bool}, s1::IntegerTypes, v2::Vec{N,T}) where {N,T<:IntegerTypes}
#     vifelse(c, vbroadcast(Vec{N,T}, s1), v2)
# end
# @inline function vifelse(c::Vec{N,Bool}, v1::Vec{N,T}, s2::IntegerTypes) where {N,T<:IntegerTypes}
#     vifelse(c, v1, vbroadcast(Vec{N,T}, s2))
# end

# for op ∈ (:muladd,)
#     rename = VECTOR_SYMBOLS[op]
#     @eval begin
#         @inline $rename(s1::IntegerTypes, v2::Vec{N,T},
#                 v3::Vec{N,T}) where {N,T<:IntegerTypes} =
#             $rename(vbroadcast(Vec{N,T}, s1), v2, v3)
#         @inline $rename(v1::Vec{N,T}, s2::IntegerTypes,
#                 v3::Vec{N,T}) where {N,T<:IntegerTypes} =
#             $rename(v1, vbroadcast(Vec{N,T}, s2), v3)
#         @inline $rename(s1::IntegerTypes, s2::IntegerTypes,
#                 v3::Vec{N,T}) where {N,T<:IntegerTypes} =
#             $rename(vbroadcast(Vec{N,T}, s1), vbroadcast(Vec{N,T}, s2), v3)
#         @inline $rename(v1::Vec{N,T}, v2::Vec{N,T},
#                 s3::IntegerTypes) where {N,T<:IntegerTypes} =
#             $rename(v1, v2, vbroadcast(Vec{N,T}, s3))
#         @inline $rename(s1::IntegerTypes, v2::Vec{N,T},
#                 s3::IntegerTypes) where {N,T<:IntegerTypes} =
#             $rename(vbroadcast(Vec{N,T}, s1), v2, vbroadcast(Vec{N,T}, s3))
#         @inline $rename(v1::Vec{N,T}, s2::IntegerTypes,
#                 s3::IntegerTypes) where {N,T<:IntegerTypes} =
#             $rename(v1, vbroadcast(Vec{N,T}, s2), vbroadcast(Vec{N,T}, s3))
#     end
# end
