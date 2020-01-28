const LOG2 = log(2)
@inline vmullog2(x::T) where {T<:Number} = vmul(T(LOG2), x) 
@inline vmullog2(v::V) where {W,T,V<:AbstractSIMDVector{W,T}} = vmul(vbroadcast(V, T(LOG2)), v)
@inline vmullog2add(x::T, y) where {T<:Number} = vfmadd(T(LOG2), x, y)
@inline vmullog2add(v::V, y) where {W,T,V<:AbstractSIMDVector{W,T}} = vfmadd(vbroadcast(V, T(LOG2)), v, y)

const LOG10 = log(10)
@inline vmullog10(x::T) where {T<:Number} = vmul(T(LOG10), x) 
@inline vmullog10(v::V) where {W,T,V<:AbstractSIMDVector{W,T}} = vmul(vbroadcast(V, T(LOG10)), v)
@inline vmullog10add(x::T, y) where {T<:Number} = vfmadd(T(LOG10), x, y)
@inline vmullog10add(v::V, y) where {W,T,V<:AbstractSIMDVector{W,T}} = vfmadd(vbroadcast(V, T(LOG10)), v, y)

const INVLOG2 = 1 / log(2)
@inline vdivlog2(x::T) where {T<:Number} = vmul(T(INVLOG2), x) 
@inline vdivlog2(v::V) where {W,T,V<:AbstractSIMDVector{W,T}} = vmul(vbroadcast(V, T(INVLOG2)), v)
@inline vdivlog2add(x::T, y) where {T<:Number} = vfmadd(T(INVLOG2), x, y)
@inline vdivlog2add(v::V, y) where {W,T,V<:AbstractSIMDVector{W,T}} = vfmadd(vbroadcast(V, T(INVLOG2)), v, y)

const INVLOG10 = 1 / log(10)
@inline vdivlog10(x::T) where {T<:Number} = vmul(T(INVLOG10), x) 
@inline vdivlog10(v::V) where {W,T,V<:AbstractSIMDVector{W,T}} = vmul(vbroadcast(V, T(INVLOG10)), v)
@inline vdivlog10add(x::T, y) where {T<:Number} = vfmadd(T(INVLOG10), x, y)
@inline vdivlog10add(v::V, y) where {W,T,V<:AbstractSIMDVector{W,T}} = vfmadd(vbroadcast(V, T(INVLOG10)), v, y)

@inline vfmaddaddone(x::T) where {T} = vfmadd(x, x, VectorizationBase.vone(promote_type(T)))
@inline vfmaddaddone(x::T1, y::T2) where {T1,T2} = vfmadd(x, y, VectorizationBase.vone(promote_type(T1,T2)))


