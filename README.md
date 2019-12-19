
`SIMDPirates.jl` is a library for `SIMD` intrinsics. The code was stolen from  [SIMD.jl](https://github.com/eschnett/SIMD.jl), whose authors and maintainers deserve credit for most of the good work here. Aside from pirating code, `SIMDPirates` also provides an `@pirate` macro that lets you imagine you're commiting type piracy:
```julua
julia> @macroexpand @pirate v1 * v2 + v3 * v4
:(SIMDPirates.vmuladd(v1, v2, SIMDPirates.vmul(v3, v4)))
```
The functions `SIMDPirates.vmuladd` and `SIMDPirates.vmul` have methods defined on `NTuple{W,Core.VecElement{<:Union{Float64,Float32}}}`. By substituting base functions for methods with appropriate definitions, you can thus use base types without actual piracy. In general, the recomended approach however is to use `SVec`, a `struct`-wrapped vector which has overloads.
```julia
julia> vbroadcast(Val(8), 0.0) + 2
SVec{8,Float64}<2.0, 2.0, 2.0, 2.0, 2.0, 2.0, 2.0, 2.0>
```

Anyone is more than welcome to take any of the code or changes I've made to this library and submit them to [SIMD.jl](https://github.com/eschnett/SIMD.jl).

### Highlights
First, generating a few random vectors (if you're interested in the SIMD generation of random numbers, please see [VectorizedRNG.jl](https://github.com/chriselrod/VectorizedRNG.jl).
```julia
julia> x = SVec(ntuple(Val(8)) do w Core.VecElement(rand()) end)
SVec{8,Float64}<0.5692277987210761, 0.33665348761817304, 0.03954926738748976, 0.3213190689556804, 0.8088511245418579, 0.35544805303664107, 0.3677375589109022, 0.4651001170793463>

julia> y = SVec(ntuple(Val(8)) do w Core.VecElement(rand()) end)
SVec{8,Float64}<0.4777741139597824, 0.06049602694925049, 0.3872217501123736, 0.486269129542215, 0.7425786365371663, 0.5857635301041568, 0.3686591983067562, 0.2412057239643277>

julia> z = SVec(ntuple(Val(8)) do w Core.VecElement(rand()) end)
SVec{8,Float64}<0.7913560950516021, 0.04884861331731183, 0.7385341388346818, 0.5228028153085258, 0.21908962866195014, 0.41415395968234314, 0.2655341712486954, 0.16997469510081653>
```
Fast flags on common operators, to allow for contractions:
```julia
julia> foo(x,y,z) = x * y - z
foo (generic function with 1 method)

julia> foo(x,y) = foo(x,y,0.0)
foo (generic function with 2 methods)
```
This results in the following asm:
```asm
#julia> @code_native debuginfo=:none foo(x,y,z)
	.text
	vmovupd	(%rsi), %zmm0
	vmovupd	(%rdx), %zmm1
	vfmsub213pd	(%rcx), %zmm0, %zmm1
	vmovapd	%zmm1, (%rdi)
	movq	%rdi, %rax
	vzeroupper
	retq
	nop

#julia> @code_native debuginfo=:none foo(x,y)
	.text
	vmovupd	(%rsi), %zmm0
	vmulpd	(%rdx), %zmm0, %zmm0
	vmovapd	%zmm0, (%rdi)
	movq	%rdi, %rax
	vzeroupper
	retq
	nopl	(%rax)
```
The arithmetic of the three argument `foo` was reduced to a single `vfmsub213pd` instruction (vectorized fused multiplication and subtraction of packed double precision numbers), while the two argument version dropped the `0.0` producing a `vmulpd` (vectorized multiplication of packed doubles).

For implementing [compensated algorithms](https://github.com/JuliaMath/AccurateArithmetic.jl), `SIMDPirates` provides functions that prevent these optimizations:
```julia
julia> efoo(x,y,z) = SIMDPirates.evsub(SIMDPirates.evmul(x, y), z)
efoo (generic function with 1 method)

julia> efoo(x,y) = efoo(x, y, x - x)
efoo (generic function with 2 methods)


julia> @code_native debuginfo=:none efoo(x,y,z)
	.text
	vmovupd	(%rsi), %zmm0
	vmulpd	(%rdx), %zmm0, %zmm0
	vsubpd	(%rcx), %zmm0, %zmm0
	vmovapd	%zmm0, (%rdi)
	movq	%rdi, %rax
	vzeroupper
	retq
	nop

julia> @code_native debuginfo=:none efoo(x,y)
	.text
	vmovupd	(%rsi), %zmm0
	vmulpd	(%rdx), %zmm0, %zmm0
	vmovapd	%zmm0, (%rdi)
	movq	%rdi, %rax
	vzeroupper
	retq
	nopl	(%rax)
```
Although it still allides subtracting (`x-x`), the multiplication and subtraction haven't contracted.

Most of the more interesting additions and changes are related to memory management.
The prefered means of masking is to use bitmasks instead of `NTuple{W,Core.VecElement{Bool}}`.
This is in large part because working with bitmasks is extremely efficient on avx512 architectures.
Note that each `Bool` is a byte, an `i8` rather than an `i1`.
Note that the zero extensions and truncations to convert back and fourth between them would almost
certainly be ellided by the compiler. The advantage of the bit representation is that it can
be convenient to generate masks by means other than a vectorized comparison. For example, if you
have a loop of length `N = 141`, and are using vectors of width 8:
```julia
julia> 141 & 7
5

julia> one(UInt8) << ans - one(UInt8)
0x1f

julia> bitstring(ans)
"00011111"

julia> x = rand(5);

julia> vload(Val(8), pointer(x), 0x1f)
SVec{8,Float64}<0.9883925090112797, 0.5963333776815305, 0.39507716254066527, 0.20452877630045485, 0.11416439490499686, 0.0, 0.0, 0.0>
```
This can be an efficient means of safely calculating the remaining iterations of a loop without segfaulting by going out of bounds.

This library also provides cartesian indexing, uses [VectorizationBase.jl](https://github.com/chriselrod/VectorizationBase.jl)'s stridedpointer.
Note that because these are pointers (rather than arrays), 0-based indices seemed more appropriate.
```julia
julia> A = randn(8,8,8);

julia> using VectorizationBase

julia> vload(Val(8), stridedpointer(A), (0,2,4))
SVec{8,Float64}<0.5652101029566953, 0.3735600961400492, -0.46186341442110324, -0.023470374325385516, -0.05667600480551983, -0.5376619417499121, -0.660267667473099, 0.4155986530326794>

julia> A[:,3,5]'
1×8 LinearAlgebra.Adjoint{Float64,Array{Float64,1}}:
 0.56521  0.37356  -0.461863  -0.0234704  -0.056676  -0.537662  -0.660268  0.415599
```
The strided pointer can also handle non-unit strides i nthe first axis:
```
julia> A = randn(8,8,8);

julia> vload(Val(8), stridedpointer(@view(A[3,:,:])), (0,4))
SVec{8,Float64}<-0.6981645894507432, -0.21264670662945478, -0.46186341442110324, 1.0916967467999321, 0.0676744641262481, 1.6641946624495672, -0.36650334364272646, -0.20951071047318678>

julia> @view(A[3,:,:])[:,5]'
1×8 LinearAlgebra.Adjoint{Float64,Array{Float64,1}}:
 -0.698165  -0.212647  -0.461863  1.0917  0.0676745  1.66419  -0.366503  -0.209511

julia> 

julia> @code_native debuginfo=:none vload(Val(8), stridedpointer(@view(A[3,:,:])), (0,4))
	.text
	movq	8(%rsi), %rax
	movq	16(%rsi), %rcx
	leaq	(,%rax,8), %r8
	imulq	(%rdx), %rax
	imulq	8(%rdx), %rcx
	addq	%rax, %rcx
	shlq	$3, %rcx
	addq	(%rsi), %rcx
	vpbroadcastq	%r8, %zmm0
	vpbroadcastq	%rcx, %zmm1
	movabsq	$139856148733824, %rax  # imm = 0x7F32CC109F80
	vpmullq	(%rax), %zmm0, %zmm0
	vpaddq	%zmm0, %zmm1, %zmm0
	kxnorw	%k0, %k0, %k1
   	vgatherqpd	(,%zmm0), %zmm1 {%k1}
	vmovapd	%zmm1, (%rdi)
	movq	%rdi, %rax
	vzeroupper
	retq
	nopw	%cs:(%rax,%rax)
```
It used a `vgatherqpd` instruction to load the unevenly spaced data. Note that this is much less efficient than a `vmov(a/u)pd`, but importantly this means that
if you write a function taking arrays as arguments, and someone passes in a view where `stride(A,1) != 1`, your function should still produce the correct answer.

Something else fun:
```julia
@inline function testcore!(ptrai::Ptr{T},ptrbi::Ptr{T},ptrci::Ptr{T},::Val{L}) where {T,L}
    ptrb = ptrbi
    ptrc = ptrci
    ptra = ptrai
    SIMDPirates.lifetime_start!(ptrai, Val(L))
    V = VectorizationBase.pick_vector(T)
    W, Wshift = VectorizationBase.pick_vector_width_shift(T)
    incr = sizeof(T) << Wshift
    for _ ∈ 1:(L>>>Wshift)
        vb = vload(V, ptrb)
        vc = vload(V, ptrc)
        vstore!(ptra, vmul(vb, vc))
        ptra += incr
        ptrb += incr
        ptrc += incr
    end
    ptra = ptrai
    out = vbroadcast(V, zero(T))
    for _ ∈ 1:(L>>>Wshift)
        out = vadd(out, vload(V, ptra))
        ptra += incr
    end
    SIMDPirates.lifetime_end!(ptrai, Val(L))
    vsum(out)
end
testsplit(d) = (pd = pointer(d); testcore!(pd, pd + 256, pd + 512, Val(32)))
```
This function takes a dot product of `b` and `c` in a stupid way:
it multiplies them elementwise, storing the results in preallocated storage `a`, before
summing up the elements of `a`.

```julia
julia> A = rand(32, 3);

julia> buff = @view(A[:,1]);

julia> b = @view(A[:,2]);

julia> c = @view(A[:,3]);

julia> using Random

julia> rand!(b); rand!(c); fill!(buff, 999.9);

julia> b' * c
9.017807879336804

julia> testsplit(A)
9.017807879336804

julia> buff'
1×32 LinearAlgebra.Adjoint{Float64,SubArray{Float64,1,Array{Float64,2},Tuple{Base.Slice{Base.OneTo{Int64}},Int64},true}}:
 999.9  999.9  999.9  999.9  999.9  999.9  999.9  999.9  999.9  999.9  999.9  999.9  999.9  999.9  999.9  999.9  999.9  999.9  999.9  999.9  999.9  999.9  999.9  999.9  999.9  999.9  999.9  999.9  999.9  999.9  999.9  999.9
```
However, the contents of `buff` are unchanged! It seems that we did not actually store into it.
```llvm
;; julia> @code_llvm debuginfo=:none testsplit(A)

define double @julia_testsplit_17617(%jl_value_t addrspace(10)* nonnull align 16 dereferenceable(40)) {
top:
  %1 = addrspacecast %jl_value_t addrspace(10)* %0 to %jl_value_t addrspace(11)*
  %2 = addrspacecast %jl_value_t addrspace(11)* %1 to %jl_value_t*
  %3 = bitcast %jl_value_t* %2 to i8**
  %4 = load i8*, i8** %3, align 8
  %5 = getelementptr i8, i8* %4, i64 256
  %6 = getelementptr i8, i8* %4, i64 512
  call void @llvm.lifetime.start.p0i8(i64 256, i8* nonnull %4)
  %ptr.i1825 = bitcast i8* %5 to <8 x double>*
  %res.i1926 = load <8 x double>, <8 x double>* %ptr.i1825, align 8
  %ptr.i1627 = bitcast i8* %6 to <8 x double>*
  %res.i1728 = load <8 x double>, <8 x double>* %ptr.i1627, align 8
  %res.i1529 = fmul reassoc nnan ninf nsz arcp contract <8 x double> %res.i1926, %res.i1728
  %7 = getelementptr i8, i8* %4, i64 576
  %8 = getelementptr i8, i8* %4, i64 320
  %ptr.i18 = bitcast i8* %8 to <8 x double>*
  %res.i19 = load <8 x double>, <8 x double>* %ptr.i18, align 8
  %ptr.i16 = bitcast i8* %7 to <8 x double>*
  %res.i17 = load <8 x double>, <8 x double>* %ptr.i16, align 8
  %res.i15 = fmul reassoc nnan ninf nsz arcp contract <8 x double> %res.i19, %res.i17
  %9 = getelementptr i8, i8* %4, i64 640
  %10 = getelementptr i8, i8* %4, i64 384
  %ptr.i18.1 = bitcast i8* %10 to <8 x double>*
  %res.i19.1 = load <8 x double>, <8 x double>* %ptr.i18.1, align 8
  %ptr.i16.1 = bitcast i8* %9 to <8 x double>*
  %res.i17.1 = load <8 x double>, <8 x double>* %ptr.i16.1, align 8
  %res.i15.1 = fmul reassoc nnan ninf nsz arcp contract <8 x double> %res.i19.1, %res.i17.1
  %11 = getelementptr i8, i8* %4, i64 704
  %12 = getelementptr i8, i8* %4, i64 448
  %ptr.i18.2 = bitcast i8* %12 to <8 x double>*
  %res.i19.2 = load <8 x double>, <8 x double>* %ptr.i18.2, align 8
  %ptr.i16.2 = bitcast i8* %11 to <8 x double>*
  %res.i17.2 = load <8 x double>, <8 x double>* %ptr.i16.2, align 8
  %res.i15.2 = fmul reassoc nnan ninf nsz arcp contract <8 x double> %res.i19.2, %res.i17.2
  %res.i12 = fadd reassoc nnan ninf nsz arcp contract <8 x double> %res.i1529, %res.i15
  %res.i12.1 = fadd reassoc nnan ninf nsz arcp contract <8 x double> %res.i12, %res.i15.1
  %res.i12.2 = fadd reassoc nnan ninf nsz arcp contract <8 x double> %res.i12.1, %res.i15.2
  call void @llvm.lifetime.end.p0i8(i64 256, i8* nonnull %4)
  %vec_4_1.i = shufflevector <8 x double> %res.i12.2, <8 x double> undef, <4 x i32> <i32 0, i32 1, i32 2, i32 3>
  %vec_4_2.i = shufflevector <8 x double> %res.i12.2, <8 x double> undef, <4 x i32> <i32 4, i32 5, i32 6, i32 7>
  %vec_4.i = fadd <4 x double> %vec_4_1.i, %vec_4_2.i
  %vec_2_1.i = shufflevector <4 x double> %vec_4.i, <4 x double> undef, <2 x i32> <i32 0, i32 1>
  %vec_2_2.i = shufflevector <4 x double> %vec_4.i, <4 x double> undef, <2 x i32> <i32 2, i32 3>
  %vec_2.i = fadd <2 x double> %vec_2_1.i, %vec_2_2.i
  %vec_1_1.i = shufflevector <2 x double> %vec_2.i, <2 x double> undef, <1 x i32> zeroinitializer
  %vec_1_2.i = shufflevector <2 x double> %vec_2.i, <2 x double> undef, <1 x i32> <i32 1>
  %vec_1.i = fadd <1 x double> %vec_1_1.i, %vec_1_2.i
  %res.i = extractelement <1 x double> %vec_1.i, i32 0
  ret double %res.i
}
```
Indeed, there are no stores. Because of the `lifetime.end`, we declared that the contents are
undefined once the function expires, therefore it does not have to write.

Unfortunately, this optimization is extremely brittle / hard to take advantage of. If there
is any possibility of aliasing, for example, it will not trigger (here, by calculating constant
offsets from a base pointer, LLVM could figure out that there was no aliasing).

