
# Type-dependent optimization flags
fastflags(::Type{T}) where {T<:IntTypes}= "nsw"
fastflags(::Type{T}) where {T<:UIntTypes} = "nuw"
fastflags(::Type{T}) where {T<:FloatingTypes} = "fast"


suffix(N::Integer, ::Type{T}) where {T<:IntegerTypes} = "v$(N)i$(8*sizeof(T))"
suffix(N::Integer, ::Type{T}) where {T<:FloatingTypes} = "v$(N)f$(8*sizeof(T))"

# Type-dependent LLVM constants
function llvmconst(::Type{T}, val) where T
    T(val) === T(0) && return "zeroinitializer"
    typ = llvmtype(T)
    "$typ $val"
end
function llvmconst(::Type{Bool}, val)
    Bool(val) === false && return "zeroinitializer"
    typ = "i1"
    "$typ $(Int(val))"
end
function llvmconst(N::Integer, ::Type{T}, val) where T
    T(val) === T(0) && return "zeroinitializer"
    typ = llvmtype(T)
    "<" * join(["$typ $val" for i in 1:N], ", ") * ">"
end
function llvmconst(N::Integer, ::Type{Bool}, val)
    Bool(val) === false && return "zeroinitializer"
    typ = "i1"
    "<" * join(["$typ $(Int(val))" for i in 1:N], ", ") * ">"
end
function llvmtypedconst(::Type{T}, val) where T
    typ = llvmtype(T)
    T(val) === T(0) && return "$typ zeroinitializer"
    "$typ $val"
end
function llvmtypedconst(::Type{Bool}, val)
    typ = "i1"
    Bool(val) === false && return "$typ zeroinitializer"
    "$typ $(Int(val))"
end

# Type-dependent LLVM intrinsics
llvmins(::Type{Val{:+}}, N, ::Type{T}) where {T <: IntegerTypes} = "add"
llvmins(::Type{Val{:-}}, N, ::Type{T}) where {T <: IntegerTypes} = "sub"
llvmins(::Type{Val{:*}}, N, ::Type{T}) where {T <: IntegerTypes} = "mul"
llvmins(::Type{Val{:div}}, N, ::Type{T}) where {T <: IntTypes} = "sdiv"
llvmins(::Type{Val{:rem}}, N, ::Type{T}) where {T <: IntTypes} = "srem"
llvmins(::Type{Val{:div}}, N, ::Type{T}) where {T <: UIntTypes} = "udiv"
llvmins(::Type{Val{:rem}}, N, ::Type{T}) where {T <: UIntTypes} = "urem"

llvmins(::Type{Val{:~}}, N, ::Type{T}) where {T <: IntegerTypes} = "xor"
llvmins(::Type{Val{:&}}, N, ::Type{T}) where {T <: IntegerTypes} = "and"
llvmins(::Type{Val{:|}}, N, ::Type{T}) where {T <: IntegerTypes} = "or"
llvmins(::Type{Val{:⊻}}, N, ::Type{T}) where {T <: IntegerTypes} = "xor"

llvmins(::Type{Val{:<<}}, N, ::Type{T}) where {T <: IntegerTypes} = "shl"
llvmins(::Type{Val{:>>>}}, N, ::Type{T}) where {T <: IntegerTypes} = "lshr"
llvmins(::Type{Val{:>>}}, N, ::Type{T}) where {T <: UIntTypes} = "lshr"
llvmins(::Type{Val{:>>}}, N, ::Type{T}) where {T <: IntTypes} = "ashr"

llvmins(::Type{Val{:(==)}}, N, ::Type{T}) where {T <: IntegerTypes} = "icmp eq"
llvmins(::Type{Val{:(!=)}}, N, ::Type{T}) where {T <: IntegerTypes} = "icmp ne"
llvmins(::Type{Val{:(>)}}, N, ::Type{T}) where {T <: IntTypes} = "icmp sgt"
llvmins(::Type{Val{:(>=)}}, N, ::Type{T}) where {T <: IntTypes} = "icmp sge"
llvmins(::Type{Val{:(<)}}, N, ::Type{T}) where {T <: IntTypes} = "icmp slt"
llvmins(::Type{Val{:(<=)}}, N, ::Type{T}) where {T <: IntTypes} = "icmp sle"
llvmins(::Type{Val{:(>)}}, N, ::Type{T}) where {T <: UIntTypes} = "icmp ugt"
llvmins(::Type{Val{:(>=)}}, N, ::Type{T}) where {T <: UIntTypes} = "icmp uge"
llvmins(::Type{Val{:(<)}}, N, ::Type{T}) where {T <: UIntTypes} = "icmp ult"
llvmins(::Type{Val{:(<=)}}, N, ::Type{T}) where {T <: UIntTypes} = "icmp ule"

llvmins(::Type{Val{:vifelse}}, N, ::Type{T}) where {T} = "select"

llvmins(::Type{Val{:+}}, N, ::Type{T}) where {T <: FloatingTypes} = "fadd"
llvmins(::Type{Val{:-}}, N, ::Type{T}) where {T <: FloatingTypes} = "fsub"
llvmins(::Type{Val{:*}}, N, ::Type{T}) where {T <: FloatingTypes} = "fmul"
llvmins(::Type{Val{:/}}, N, ::Type{T}) where {T <: FloatingTypes} = "fdiv"
llvmins(::Type{Val{:inv}}, N, ::Type{T}) where {T <: FloatingTypes} = "fdiv"
llvmins(::Type{Val{:rem}}, N, ::Type{T}) where {T <: FloatingTypes} = "frem"

llvmins(::Type{Val{:(==)}}, N, ::Type{T}) where {T <: FloatingTypes} = "fcmp oeq"
llvmins(::Type{Val{:(!=)}}, N, ::Type{T}) where {T <: FloatingTypes} = "fcmp une"
llvmins(::Type{Val{:(>)}}, N, ::Type{T}) where {T <: FloatingTypes} = "fcmp ogt"
llvmins(::Type{Val{:(>=)}}, N, ::Type{T}) where {T <: FloatingTypes} = "fcmp oge"
llvmins(::Type{Val{:(<)}}, N, ::Type{T}) where {T <: FloatingTypes} = "fcmp olt"
llvmins(::Type{Val{:(<=)}}, N, ::Type{T}) where {T <: FloatingTypes} = "fcmp ole"

llvmins(::Type{Val{:^}}, N, ::Type{T}) where {T <: FloatingTypes} =
    "@llvm.pow.$(suffix(N,T))"
llvmins(::Type{Val{:abs}}, N, ::Type{T}) where {T<:FloatingTypes} =
    "@llvm.fabs.$(suffix(N,T))"
llvmins(::Type{Val{:ceil}}, N, ::Type{T}) where {T<:FloatingTypes} =
    "@llvm.ceil.$(suffix(N,T))"
llvmins(::Type{Val{:copysign}}, N, ::Type{T}) where {T<:FloatingTypes} =
    "@llvm.copysign.$(suffix(N,T))"
llvmins(::Type{Val{:cos}}, N, ::Type{T}) where {T<:FloatingTypes} =
    "@llvm.cos.$(suffix(N,T))"
llvmins(::Type{Val{:exp}}, N, ::Type{T}) where {T<:FloatingTypes} =
    "@llvm.exp.$(suffix(N,T))"
llvmins(::Type{Val{:exp2}}, N, ::Type{T}) where {T<:FloatingTypes} =
    "@llvm.exp2.$(suffix(N,T))"
llvmins(::Type{Val{:floor}}, N, ::Type{T}) where {T<:FloatingTypes} =
    "@llvm.floor.$(suffix(N,T))"
llvmins(::Type{Val{:fma}}, N, ::Type{T}) where {T<:FloatingTypes} =
    "@llvm.fma.$(suffix(N,T))"
llvmins(::Type{Val{:log}}, N, ::Type{T}) where {T<:FloatingTypes} =
    "@llvm.log.$(suffix(N,T))"
llvmins(::Type{Val{:log10}}, N, ::Type{T}) where {T<:FloatingTypes} =
    "@llvm.log10.$(suffix(N,T))"
llvmins(::Type{Val{:log2}}, N, ::Type{T}) where {T<:FloatingTypes} =
    "@llvm.log2.$(suffix(N,T))"
llvmins(::Type{Val{:max}}, N, ::Type{T}) where {T<:FloatingTypes} =
    "@llvm.maxnum.$(suffix(N,T))"
llvmins(::Type{Val{:min}}, N, ::Type{T}) where {T<:FloatingTypes} =
    "@llvm.minnum.$(suffix(N,T))"
llvmins(::Type{Val{:muladd}}, N, ::Type{T}) where {T<:FloatingTypes} =
    "@llvm.fmuladd.$(suffix(N,T))"
llvmins(::Type{Val{:powi}}, N, ::Type{T}) where {T<:FloatingTypes} =
    "@llvm.powi.$(suffix(N,T))"
llvmins(::Type{Val{:round}}, N, ::Type{T}) where {T<:FloatingTypes} =
    "@llvm.rint.$(suffix(N,T))"
llvmins(::Type{Val{:sin}}, N, ::Type{T}) where {T<:FloatingTypes} =
    "@llvm.sin.$(suffix(N,T))"
llvmins(::Type{Val{:sqrt}}, N, ::Type{T}) where {T<:FloatingTypes} =
    "@llvm.sqrt.$(suffix(N,T))"
llvmins(::Type{Val{:trunc}}, N, ::Type{T}) where {T<:FloatingTypes} =
    "@llvm.trunc.$(suffix(N,T))"

# Convert between LLVM scalars, vectors, and arrays

function scalar2vector(vec, siz, typ, sca)
    instrs = []
    accum(nam, i) = i<0 ? "undef" : i==siz-1 ? nam : "$(nam)_iter$i"
    for i in 0:siz-1
        push!(instrs,
            "$(accum(vec,i)) = " *
                "insertelement <$siz x $typ> $(accum(vec,i-1)), " *
                "$typ $sca, i32 $i")
    end
    instrs
end

function array2vector(vec, siz, typ, arr, tmp="$(arr)_av")
    instrs = []
    accum(nam, i) = i<0 ? "undef" : i==siz-1 ? nam : "$(nam)_iter$i"
    for i in 0:siz-1
        push!(instrs, "$(tmp)_elem$i = extractvalue [$siz x $typ] $arr, $i")
        push!(instrs,
            "$(accum(vec,i)) = " *
                "insertelement <$siz x $typ> $(accum(vec,i-1)), " *
                "$typ $(tmp)_elem$i, i32 $i")
    end
    instrs
end

function vector2array(arr, siz, typ, vec, tmp="$(vec)_va")
    instrs = []
    accum(nam, i) = i<0 ? "undef" : i==siz-1 ? nam : "$(nam)_iter$i"
    for i in 0:siz-1
        push!(instrs,
            "$(tmp)_elem$i = extractelement <$siz x $typ> $vec, i32 $i")
        push!(instrs,
            "$(accum(arr,i)) = "*
                "insertvalue [$siz x $typ] $(accum(arr,i-1)), " *
                "$typ $(tmp)_elem$i, $i")
    end
    instrs
end

# TODO: change argument order
function subvector(vec, siz, typ, rvec, rsiz, roff, tmp="$(rvec)_sv")
    instrs = []
    accum(nam, i) = i<0 ? "undef" : i==rsiz-1 ? nam : "$(nam)_iter$i"
    @assert 0 <= roff
    @assert roff + rsiz <= siz
    for i in 0:rsiz-1
        push!(instrs,
            "$(tmp)_elem$i = extractelement <$siz x $typ> $vec, i32 $(roff+i)")
        push!(instrs,
            "$(accum(rvec,i)) = " *
                "insertelement <$rsiz x $typ> $(accum(rvec,i-1)), " *
                "$typ $(tmp)_elem$i, i32 $i")
    end
    instrs
end

function extendvector(vec, siz, typ, voff, vsiz, val, rvec, tmp="$(rvec)_ev")
    instrs = []
    accum(nam, i) = i<0 ? "undef" : i==siz+vsiz-1 ? nam : "$(nam)_iter$i"
    rsiz = siz + vsiz
    for i in 0:siz-1
        push!(instrs,
            "$(tmp)_elem$i = extractelement <$siz x $typ> $vec, i32 $i")
        push!(instrs,
            "$(accum(rvec,i)) = " *
                "insertelement <$rsiz x $typ> $(accum(rvec,i-1)), " *
                "$typ $(tmp)_elem$i, i32 $i")
    end
    for i in siz:siz+vsiz-1
        push!(instrs,
            "$(accum(rvec,i)) = " *
                "insertelement <$rsiz x $typ> $(accum(rvec,i-1)), $val, i32 $i")
    end
    instrs
end

# Element-wise access

# export setindex
@generated function setindex(v::Vec{N,T}, x::Number, ::Type{Val{I}}) where {N,T,I}
    @assert isa(I, Integer)
    1 <= I <= N || throw(BoundsError())
    typ = llvmtype(T)
    ityp = llvmtype(Int)
    vtyp = "<$N x $typ>"
    decls = []
    instrs = []
    push!(instrs, "%res = insertelement $vtyp %0, $typ %1, $ityp $(I-1)")
    push!(instrs, "ret $vtyp %res")
    quote
        $(Expr(:meta, :inline))
        Base.llvmcall($((join(decls, "\n"), join(instrs, "\n"))),
            NTuple{N,VE{T}}, Tuple{NTuple{N,VE{T}}, T}, v.elts, T(x))
    end
end

@generated function setindex(v::Vec{N,T}, x::Number, i::Int) where {N,T}
    typ = llvmtype(T)
    ityp = llvmtype(Int)
    vtyp = "<$N x $typ>"
    decls = []
    instrs = []
    push!(instrs, "%res = insertelement $vtyp %0, $typ %2, $ityp %1")
    push!(instrs, "ret $vtyp %res")
    quote
        $(Expr(:meta, :inline))
        @boundscheck 1 <= i <= N || throw(BoundsError())
        Base.llvmcall($((join(decls, "\n"), join(instrs, "\n"))),
            NTuple{N,VE{T}}, Tuple{NTuple{N,VE{T}}, Int, T},
            v.elts, i-1, T(x))
    end
end
setindex(v::Vec{N,T}, x::Number, i) where {N,T} = setindex(v, x, Int(i))
# Type conversion

@generated function pirate_reinterpret(::Type{Vec{N,R}},
        v1::Vec{N1,T1}) where {N,R,N1,T1}
    if N*sizeof(R) != N1*sizeof(T1)
        throw("N*sizeof(R) == N1*sizeof(T1) is not true; Trying to reinterpret to a size of $N * $(sizeof(R)) from a size of $N1 * $(sizeof(T1))")
    end
    # @assert N*sizeof(R) == N1*sizeof(T1)
    typ1 = llvmtype(T1)
    vtyp1 = "<$N1 x $typ1>"
    typr = llvmtype(R)
    vtypr = "<$N x $typr>"
    decls = String[]
    instrs = String[]
    push!(instrs, "%res = bitcast $vtyp1 %0 to $vtypr")
    push!(instrs, "ret $vtypr %res")
    quote
        $(Expr(:meta, :inline))
        Base.llvmcall($((join(decls, "\n"), join(instrs, "\n"))),
            Vec{$N,$R}, Tuple{Vec{$N1,$T1}}, v1)
    end
end
@generated function pirate_reinterpret(::Type{Vec{N,R}},
        v1::Vec{N1,UInt128}) where {N,R,N1}
    T1 = UInt128
    @assert N*sizeof(R) == N1*sizeof(T1)
    typ1 = llvmtype(T1)
    vtyp1 = "[$N1 x $typ1]"
    typr = llvmtype(R)
    vtypr = "<$N x $typr>"
    decls = String[]
    instrs = String[]
    push!(instrs, "%res = bitcast $vtyp1 %0 to $vtypr")
    push!(instrs, "ret $vtypr %res")
    quote
        $(Expr(:meta, :inline))
        Base.llvmcall($((join(decls, "\n"), join(instrs, "\n"))),
            Vec{$N,$R}, Tuple{Vec{$N1,$T1}}, v1)
    end
end
@inline function Base.reinterpret(::Type{SVec{N,R}}, v1::AbstractStructVec{N1,T1}) where {N,R,N1,T1}
    SVec(pirate_reinterpret(Vec{N,R}, extract_data(v1)))
end
@inline function Base.reinterpret(::Type{SVec{N,R}}, v1::AbstractStructVec{N1,UInt128}) where {N,R,N1}
    SVec(pirate_reinterpret(Vec{N,R}, extract_data(v1)))
end

const FASTOPS = Set((:+, :-, :*, :/, :log, :log2, :log10, :exp, :exp2, :exp10, :muladd, :fma))

const VECTOR_SYMBOLS = Dict{Symbol,Symbol}(
    :(==) => :visequal,
    :(!=) => :vnot_equal,
    :(<) => :vless,
    :(<=) => :vless_or_equal,
    :(>) => :vgreater,
    :(>=) => :vgreater_or_equal,
    :(<<) => :vleft_bitshift,
    :(>>) => :vright_bitshift,
    :(>>>) => :vuright_bitshift,
    :(&) => :vand,
    :(|) => :vor,
    :(⊻) => :vxor,
    :(+) => :vadd,
    :(-) => :vsub,
    :(*) => :vmul,
    :(/) => :vfdiv,
    :(÷) => :vidiv,
    :(%) => :vrem,
    :div => :vdiv,
    :rem => :vrem,
    :(~) => :vbitwise_not,
    :(!) => :vnot,
    :(^) => :vpow,
    :abs => :vabs,
    :floor => :vfloor,
    :ceil => :vceil,
    :round => :vround,
    :sin => :vsin,
    :cos => :vcos,
    :exp => :vexp,
    :exp2 => :vexp2,
    :exp10 => :vexp10,
    :inv => :vinv,
    :log => :vlog,
    :log10 => :vlog10,
    :log2 => :vlog2,
    :sqrt => :vsqrt,
    :trunc => :vtrunc,
    :sign => :vsign,
    :copysign => :vcopysign,
    :flipsign => :vflipsign,
    :max => :vmax,
    :min => :vmin,
    :fma => :vfma,
    :muladd => :vmuladd,
    :all => :vall,
    :any => :vany,
    :maximum => :vmaximum,
    :minimum => :vminimum,
    :prod => :vprod,
    :sum => :vsum,
    :reduce => :vreduce,
    :isfinite => :visfinite,
    :isinf => :visinf,
    :isnan => :visnan,
    :issubnormal => :vissubnormal,
    :fmadd => :vfmadd,
    :fmsub => :vfmsub,
    :fnmadd => :vfnmadd,
    :fnmsub => :vfnmsub
)
