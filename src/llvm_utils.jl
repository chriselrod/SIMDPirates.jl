
# Type-dependent optimization flags
function fastflags(@nospecialize(T))
    if T <: IntTypes
        s = "nsw"
    elseif T <: UIntTypes
        s = "nuw"
    else#if T <: FloatingTypes
        s = "fast"
        # s = "nnan ninf nsz arcp contract reassoc"
        # s = "nnan ninf nsz arcp contract reassoc"
    end
    return s
end

function suffix(N::Int, @nospecialize(T))
    if T <: Ptr
        if T <: Ptr{<:IntegerTypes}
            t = "p0i"
        else#if T <: Ptr{<:FloatingTypes}
            t = "p0f"
        end
    elseif T <: IntegerTypes
        t = "i"
    else#if T <: FloatingTypes
        t = "f"
    end
    "v$(N)$(t)$(8sizeof(T))"
end


# Type-dependent LLVM constants
function llvmconst(T, val)
    iszero(val) && return "zeroinitializer"
    typ = llvmtype(T)
    "$typ $val"
end
function llvmconst(::Type{Bool}, val)
    Bool(val) || return "zeroinitializer"
    typ = "i1"
    "$typ $(Int(val))"
end
function llvmconst(N::Integer, T, val)
    isa(val, Number) && iszero(val) && return "zeroinitializer"
    typ = llvmtype(T)
    "<" * join(["$typ $val" for i in 1:N], ", ") * ">"
end
function llvmconst(N::Integer, ::Type{Bool}, val)
    Bool(val) || return "zeroinitializer"
    typ = "i1"
    "<" * join(["$typ $(Int(val))" for i in 1:N], ", ") * ">"
end
function llvmtypedconst(T, val)
    typ = llvmtype(T)
    iszero(val) && return "$typ zeroinitializer"
    "$typ $val"
end
function llvmtypedconst(::Type{Bool}, val)
    typ = "i1"
    Bool(val) || return "$typ zeroinitializer"
    "$typ $(Int(val))"
end

# Type-dependent LLVM intrinsics
const LLVM_INS_Integer = Dict{Symbol,String}()
const LLVM_INS_Int = Dict{Symbol,String}()
const LLVM_INS_UInt = Dict{Symbol,String}()
const LLVM_INS = Dict{Symbol,String}()
const LLVM_OVERLOADED_INS = Dict{Symbol,String}()



LLVM_INS_Integer[:+] = "add"
LLVM_INS_Integer[:-] = "sub"
LLVM_INS_Integer[:*] = "mul"
LLVM_INS_Int[:÷] = "sdiv"
LLVM_INS_Int[:div] = "sdiv"
LLVM_INS_Int[:rem] = "srem"
LLVM_INS_Int[:%] = "srem"
LLVM_INS_UInt[:÷] = "udiv"
LLVM_INS_UInt[:div] = "udiv"
LLVM_INS_UInt[:rem] = "urem"
LLVM_INS_UInt[:%] = "urem"

LLVM_INS_Integer[:~] = "xor"
LLVM_INS_Integer[:&] = "and"
LLVM_INS_Integer[:|] = "or"
LLVM_INS_Integer[:⊻] = "xor"

LLVM_INS_Integer[:<<] = "shl"
LLVM_INS_Integer[:>>>] = "lshr"
LLVM_INS_UInt[:>>] = "lshr"
LLVM_INS_Int[:>>] = "ashr"

LLVM_INS_Integer[:(==)] = "icmp eq"
LLVM_INS_Integer[:(!=)] = "icmp ne"
LLVM_INS_Int[:(>)] = "icmp sgt"
LLVM_INS_Int[:(>=)] = "icmp sge"
LLVM_INS_Int[:(<)] = "icmp slt"
LLVM_INS_Int[:(<=)] = "icmp sle"
LLVM_INS_UInt[:(>)] = "icmp ugt"
LLVM_INS_UInt[:(>=)] = "icmp uge"
LLVM_INS_UInt[:(<)] = "icmp ult"
LLVM_INS_UInt[:(<=)] = "icmp ule"

LLVM_INS[:vifelse] = "select"

LLVM_INS[:+] = "fadd"
LLVM_INS[:-] = "fsub"
LLVM_INS[:*] = "fmul"
LLVM_INS[:/] = "fdiv"
LLVM_INS[:inv] = "fdiv"
LLVM_INS[:rem] = "frem"

LLVM_INS[:(==)] = "fcmp contract reassoc nsz oeq"
LLVM_INS[:(!=)] = "fcmp fast une"
LLVM_INS[:(>)] = "fcmp fast ogt"
LLVM_INS[:(>=)] = "fcmp fast oge"
LLVM_INS[:(<)] = "fcmp fast olt"
LLVM_INS[:(<=)] = "fcmp fast ole"

LLVM_OVERLOADED_INS[:^] = "@llvm.pow."
LLVM_OVERLOADED_INS[:abs] = "@llvm.fabs."
LLVM_OVERLOADED_INS[:ceil] = "@llvm.ceil."
LLVM_OVERLOADED_INS[:copysign] = "@llvm.copysign."
LLVM_OVERLOADED_INS[:cos] = "@llvm.cos."
LLVM_OVERLOADED_INS[:exp] = "@llvm.exp."
LLVM_OVERLOADED_INS[:exp2] = "@llvm.exp2."
LLVM_OVERLOADED_INS[:floor] = "@llvm.floor."
LLVM_OVERLOADED_INS[:fma] = "@llvm.fma."
LLVM_OVERLOADED_INS[:log] = "@llvm.log."
LLVM_OVERLOADED_INS[:log10] = "@llvm.log10."
LLVM_OVERLOADED_INS[:log2] = "@llvm.log2."
LLVM_OVERLOADED_INS[:max] = "@llvm.maxnum."
LLVM_OVERLOADED_INS[:min] = "@llvm.minnum."
LLVM_OVERLOADED_INS[:muladd] = "@llvm.fmuladd."
LLVM_OVERLOADED_INS[:powi] = "@llvm.powi."
LLVM_OVERLOADED_INS[:round] = "@llvm.rint."
LLVM_OVERLOADED_INS[:sin] = "@llvm.sin."
LLVM_OVERLOADED_INS[:sqrt] =  "@llvm.sqrt."
LLVM_OVERLOADED_INS[:trunc] = "@llvm.trunc."

function llvmins(func::Symbol, N::Int, T)::String
    ins = get(LLVM_OVERLOADED_INS, func, nothing)
    ins === nothing || return ins * suffix(N, T)
    if T <: IntegerTypes
        d = T <: IntTypes ? LLVM_INS_Int : LLVM_INS_UInt
        return get(d, func) do
            LLVM_INS_Integer[func]
        end
    else
        return LLVM_INS[func]
    end
end


# Convert between LLVM scalars, vectors, and arrays

function scalar2vector(vec, siz, typ, sca)
    tempvec = vec * "_temporary"
    String[
        "$tempvec = insertelement <$siz x $typ> undef, $typ $sca, i32 0",
        "$vec = shufflevector <$siz x $typ> $tempvec, <$siz x $typ> undef, <$siz x i32> zeroinitializer"
    ]
end

function array2vector(vec, siz, typ, arr, tmp="$(arr)_av")
    instrs = String[]
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
    instrs = String[]
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
    instrs = String[]
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
    instrs = String[]
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

# Type conversion

@generated function vreinterpret(::Type{_Vec{_W,R}}, v1::_Vec{_W1,T1}) where {_W,R,_W1,T1}
    W = _W + 1; W1 = _W1 + 1
    if W*sizeof(R) != W1*sizeof(T1)
        throw("W*sizeof(R) == W1*sizeof(T1) is not true; Trying to reinterpret to a size of $W * $(sizeof(R)) from a size of $W1 * $(sizeof(T1))")
    end
    # @assert W*sizeof(R) == W1*sizeof(T1)
    typ1 = llvmtype(T1)
    vtyp1 = "<$W1 x $typ1>"
    typr = R <: Ptr ? llvmtype(Int) : llvmtype(R)
    vtypr = "<$W x $typr>"
    decls = String[]
    instrs = String[]
    push!(instrs, "%res = bitcast $vtyp1 %0 to $vtypr")
    push!(instrs, "ret $vtypr %res")
    quote
        $(Expr(:meta, :inline))
        llvmcall($((join(decls, "\n"), join(instrs, "\n"))),
            Vec{$W,$R}, Tuple{Vec{$W1,$T1}}, v1)
    end
end
@generated function vreinterpret(::Type{_Vec{_W,R}}, v1::_Vec{_W1,UInt128}) where {_W,R,_W1}
    W = _W + 1; W1 = _W1 + 1
    T1 = UInt128
    @assert W*sizeof(R) == W1*sizeof(T1)
    typ1 = llvmtype(T1)
    vtyp1 = "[$W1 x $typ1]"
    typr = R <: Ptr ? llvmtype(Int) : llvmtype(R)
    vtypr = "<$W x $typr>"
    decls = String[]
    instrs = String[]
    push!(instrs, "%res = bitcast $vtyp1 %0 to $vtypr")
    push!(instrs, "ret $vtypr %res")
    quote
        $(Expr(:meta, :inline))
        llvmcall($((join(decls, "\n"), join(instrs, "\n"))),
            Vec{$W,$R}, Tuple{Vec{$W1,$T1}}, v1)
    end
end
@inline function Base.reinterpret(::Type{SVec{W,R}}, v1::SVec{W1,T1}) where {W,R,W1,T1}
    SVec(vreinterpret(Vec{W,R}, extract_data(v1)))
end
@inline function Base.reinterpret(::Type{SVec{W,R}}, v1::SVec{W1,UInt128}) where {W,R,W1}
    SVec(vreinterpret(Vec{W,R}, extract_data(v1)))
end
@inline vrem(v::_Vec{W,I}, ::Type{_Vec{W,U}}) where {W,I<:Signed,U<:Unsigned} = vrem(vreinterpret(_Vec{W,U}, v), _Vec{W,U})
@inline vrem(v::_Vec{W,U}, ::Type{_Vec{W,I}}) where {W,I<:Signed,U<:Unsigned} = vrem(vreinterpret(_Vec{W,I}, v), _Vec{W,I})

@inline function assume(b::Bool)
    llvmcall(("declare void @llvm.assume(i1)", """
    %b = trunc i8 %0 to i1
    call void @llvm.assume(i1 %b)
    ret void
    """), Nothing, Tuple{Bool}, b)
end
@inline function expect(b::Bool)
    llvmcall(("declare i1 @llvm.expect.i1(i1, i1)", """
    %b = trunc i8 %0 to i1
    %actual = call i1 @llvm.expect.i1(i1 %b, i1 true)
    %byte = zext i1 %actual to i8
    ret i8 %byte
    """), Bool, Tuple{Bool}, b)
end
@generated function expect(i::I, ::Val{N}) where {I <: Integer, N}
    ityp = 'i' * string(8sizeof(I))
    decls = "declare $ityp @llvm.expect.$ityp($ityp, $ityp)"
    instrs = """
    %actual = call $ityp @llvm.expect.$ityp($ityp %0, $ityp $N)
    ret $ityp %actual
    """
    quote
        $(Expr(:meta,:inline))
        llvmcall($((decls, instrs)), $I, Tuple{$I}, i)
    end
end

const FASTOPS = Set((:+, :-, :*, :/, :log, :log2, :log10, :exp, :exp2, :exp10, :sqrt, :pow, :sin, :cos))#, :inv, :muladd, :fma

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
    :(÷) => :vdiv,
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
