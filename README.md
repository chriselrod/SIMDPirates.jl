
`SIMDPirates.jl` is a library for `SIMD` intrinsics. The code was stolen from  [SIMD.jl](https://github.com/eschnett/SIMD.jl), whose authors and maintainers deserve credit for the good work here. Aside from pirating code, `SIMDPirates` threatens to pirate methods. `SIMD.jl` provides methods for a custom vector type defined as:
```julia
const VE = Base.VecElement
struct Vec{N,T<:ScalarTypes} <: DenseArray{T,1}   # <: Number
    elts::NTuple{N,VE{T}}
    @inline Vec{N,T}(elts::NTuple{N, VE{T}}) where {N,T} = new{N,T}(elts)
end
```
This allows `SIMD` to extend methods for functions such as `+`, `*`, or `sqrt` to dispatch to the appropriate `SIMD` intrinsic methods when called on the `SIMD.Vec` type. However, SIMDPirates instead defines an alias to a base type:
```julia
const Vec{N,T} = NTuple{N,Core.VecElement{T}}
```
instead of wrapping the `NTuple{N,Core.VecElement{T}}` in another `struct` type. This means that `SIMDPirates` cannot extend these base methods without committing type piracy. For this reason, base methods are not extended. Instead, the library provides unexported methods, generally prefixed with a `v`. Instead of `a + b`, we have `SIMDPirates.vadd(a,b)`. In place of `fma(a,b,c)`, it is `SIMDPirates.vfma(a,b,c)`, etc.

Type piracy is extending methods of functions defined outside of the present library to types defined outside the library. This is not kosher, because then the act of importing a library suddenly changes global state; any other code that assumes different behavior between these functions and types now suddenly has bugs. However, the `@pirate` macro allows you to locally pirate only the methods in the following expression, avoiding the dangers of piracy.
```julia
julia> using SIMDPirates

julia> C = randn(8,8);

julia> v1 = vload(Vec{8,Float64},C, 1);

julia> v2 = vload(Vec{8,Float64},C, 9);

julia> v3 = vload(Vec{8,Float64},C,17);

julia> v4 = vload(Vec{8,Float64},C,25);

julia> @macroexpand @pirate sqrt(v1)
:((SIMDPirates.vsqrt)(v1))

julia> @macroexpand @pirate v1 + v2
:((SIMDPirates.vadd)(v1, v2))

julia> @macroexpand @pirate v1 + v2 * v3
:((SIMDPirates.vmuladd)(v2, v3, v1))

julia> @macroexpand @pirate v1 * v2 + v3 * v4
:((SIMDPirates.vmuladd)(v1, v2, (SIMDPirates.vmul)(v3, v4)))
```

Why the trouble?
LLVM has an easier time optimizing code without the wrapper. Using `SIMDPirates.jl`:
```julia
julia> @code_native SIMDPirates.vmuladd(v1, v2, v3)
	.text
; Function vmuladd {
; Location: floating_point_arithmetic.jl:48
; Function llvmwrap; {
; Location: llvmwrap.jl:148
; Function llvmwrap; {
; Location: llvmwrap.jl:148
; Function macro expansion; {
; Location: floating_point_arithmetic.jl:48
	vfmadd213pd	%zmm2, %zmm1, %zmm0
;}}}
	retq
	nopw	(%rax,%rax)
;}
```
while `SIMD.jl` results in plenty of unnecessary mov instructions:
```julia
julia> @code_native muladd(v1, v2, v3)
	.text
; Function muladd {
; Location: SIMD.jl:1007
; Function llvmwrap; {
; Location: SIMD.jl:663
; Function llvmwrap; {
; Location: SIMD.jl:663
; Function macro expansion; {
; Location: SIMD.jl:1007
	vmovupd	(%rsi), %zmm0
	vmovupd	(%rdx), %zmm1
	vfmadd213pd	(%rcx), %zmm0, %zmm1
;}}}
	vmovapd	%zmm1, (%rdi)
	movq	%rdi, %rax
	vzeroupper
	retq
	nop
;}
```
Compared to `SIMDPirates.jl`, `SIMD.jl` resulted in lots of unnecessary instructions; in particular: two `vmovupd`, one `vmovapd`, and one `movq`.
The compiler can often eliminate many of these instructions, but rarely all $-$ as seen here.

Additionally, `SIMDPirates.jl` works better with `SLEEFwrap.jl` for high performance vectorized elementary functions. `SIMD.jl` results in several extra move instructions per call, and in the case of avx-512, has been highly unstable with frequent segmentation faults (although `SLEEFwrap.sincos` still segfaults in benchmarks).
