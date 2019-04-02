#ifndef DX_MATH_MERCURY_H
#define DX_MATH_MERCURY_H
#pragma once

#include "dxmath.h"

/* Dummy type to ensure we don't mix wrapped and unwrapped vectors. */
struct MerDX_VECTOR4_wrapper;
struct MerDX_MATRIX44_wrapper;

/* macros allow us to assign over an existing vector. These are used,
 * in combination with clobbering inputs, to avoid as many allocations as
 * possible.
 */

#ifdef DXMATH_NO_INTRINSICS

/* No need to bother with alignment when no intrinsics are being used. */
#define MerDX_Allocate(TYPE) MR_GC_malloc_atomic(sizeof(TYPE))
#define MerDX_ComputeOffset(X) ((MR_Word)(X))

#else

/* Ensure alignment when using intrinsics. */
#define MerDX_Allocate(TYPE) \
    MR_GC_malloc_atomic(MR_round_up(sizeof(TYPE), 16))

MR_Word DXMATH_CONSTEXPR DXMATH_CALL MerDX_ComputeOffset(const void *const w){
    MR_Word data = (MR_Word)w;
    if((data & 0xF) == 0) return data;
    data = (data-1) & (~((MR_Word)0xF));
    data += 0x10;
    return data;
}

#endif

#define MerDX_AllocateVector4() \
   ((struct MerDX_VECTOR4_wrapper*) MerDX_Allocate(DXMATH_VECTOR4))

#define MerDX_AllocateMatrix44() \
   ((struct MerDX_MATRIX44_wrapper*) MerDX_Allocate(DXMATH_MATRIX44))
    
typedef DXMATH_VECTOR4* MerDX_DXMATH_VECTOR4_PTR;
typedef DXMATH_MATRIX44* MerDX_DXMATH_MATRIX44_PTR;

MerDX_DXMATH_VECTOR4_PTR DXMATH_CONSTEXPR DXMATH_CALL
MerDX_Vector4Offset(struct MerDX_VECTOR4_wrapper *w){
    return (DXMATH_VECTOR4*)MerDX_ComputeOffset(w);
}

MerDX_DXMATH_MATRIX44_PTR DXMATH_CONSTEXPR DXMATH_CALL
MerDX_Matrix44Offset(struct MerDX_MATRIX44_wrapper *w){
    return (DXMATH_MATRIX44*)MerDX_ComputeOffset(w);
}

typedef struct MerDX_VECTOR4_wrapper* MerDX_VECTOR4_PTR;
typedef struct MerDX_MATRIX44_wrapper* MerDX_MATRIX44_PTR;

MerDX_VECTOR4_PTR DXMATH_CONSTEXPR DXMATH_CALL
MerDX_CreateVector4Assign(DXMATH_VECTOR4 V){
    struct MerDX_VECTOR4_wrapper *const wrapper = MerDX_AllocateVector4();
    *MerDX_Vector4Offset(wrapper) = V;
    return wrapper;
}

MerDX_MATRIX44_PTR DXMATH_CONSTEXPR DXMATH_CALL
MerDX_CreateMatrix44Assign(DXMATH_MATRIX44 V){
    struct MerDX_MATRIX44_wrapper *const wrapper = MerDX_AllocateMatrix44();
    *MerDX_Matrix44Offset(wrapper) = V;
    return wrapper;
}

#define MerDX_CreateVector4(FX, FY, FZ, FW) \
    MerDX_CreateVector4Assign(DXMATH_Vector4Set((FX), (FY), (FZ), (FW)))

#define MerDX_CreateMatrix44(AX, AY, AZ, AW, BX, BY, BZ, BW, CX, CY, CZ, CW, DX, DY, DZ, DW) \
    MerDX_CreateMatrix44Assign(XMMatrix44Set( \
        (AX), (AY), (AZ), (AW), \
        (BX), (BY), (BZ), (BW), \
        (CX), (CY), (CZ), (CW), \
        (DX), (DY), (DZ), (DW)))

#define MerDX_Vector4Unwrap(VEC_WRAP) \
    (*MerDX_Vector4Offset((VEC_WRAP)))

#define MerDX_Matrix44Unwrap(MAT_WRAP) \
    (*MerDX_Matrix44Offset((MAT_WRAP)))

#define MerDX_DESTRUCTIVE_BINOP(MERDX_OP, MERDX_OVERWRITE, MERDX_V1, MERDX_V2) \
    (MERDX_OVERWRITE); \
    do{ \
    const DXMATH_VECTOR4 MerDX_DESTRUCTIVE_BINOP_vec = \
        MERDX_OP(MerDX_Vector4Unwrap((MERDX_V1)), MerDX_Vector4Unwrap((MERDX_V2))); \
    MerDX_Vector4Unwrap(Out) = MerDX_DESTRUCTIVE_BINOP_vec; \
    } while(0)

#endif /* DX_MATH_MERCURY_H */
