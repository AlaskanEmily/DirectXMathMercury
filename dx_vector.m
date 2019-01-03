% Copyright (c) 2019 Alaskan Emily, Transnat Games
%
% This software is provided 'as-is', without any express or implied warranty.
% In no event will the authors be held liable for any damages arising from
% the use of this software.
%
% Permission is granted to anyone to use this software for any purpose,
% including commercial applications, and to alter it and redistribute it
% freely, subject to the following restrictions:
%
%   1. The origin of this software must not be misrepresented; you must not
%      claim that you wrote the original software. If you use this software
%      in a product, an acknowledgment in the product documentation would be
%      appreciated but is not required.
%
%   2. Altered source versions must be plainly marked as such, and must not
%      be misrepresented as being the original software.
%
%   3. This notice may not be removed or altered from any source distribution.
%

:- module dx_vector.
%=============================================================================%
% Bindings to DirectX Math vector4.
% Provides pure Mercury fallbacks for grades that do not interface with C.
:- interface.
%=============================================================================%

:- type vector. % ---> vector(x::float, y::float, z::float, w::float).

%-----------------------------------------------------------------------------%
% Used to simulate a vector functor.
:- func vector(float, float, float, float) = vector.
:- mode vector(in, in, in, in) = uo is det.
:- mode vector(uo, uo, uo, uo) = in is det.

%-----------------------------------------------------------------------------%
% Field access functions.
:- func x(vector) = float.
:- func y(vector) = float.
:- func z(vector) = float.
:- func w(vector) = float.

%-----------------------------------------------------------------------------%

:- func 'x :='(vector, float) = vector.
:- func 'y :='(vector, float) = vector.
:- func 'z :='(vector, float) = vector.
:- func 'w :='(vector, float) = vector.

%-----------------------------------------------------------------------------%
% Arithmetic operations. These are multi-moded to allow destructive update,
% limiting the number of allocations that must occur.

:- func (vector) + (vector) = (vector).
:- mode (in) + (in) = (uo) is det.
:- mode (di) + (in) = (uo) is det.
:- mode (in) + (di) = (uo) is det.
:- mode (di) + (di) = (uo) is det.

%-----------------------------------------------------------------------------%

:- func (vector) - (vector) = (vector).
:- mode (in) - (in) = (uo) is det.
:- mode (di) - (in) = (uo) is det.
:- mode (in) - (di) = (uo) is det.
:- mode (di) - (di) = (uo) is det.

%-----------------------------------------------------------------------------%

:- func (vector) * (vector) = (vector).
:- mode (in) * (in) = (uo) is det.
:- mode (di) * (in) = (uo) is det.
:- mode (in) * (di) = (uo) is det.
:- mode (di) * (di) = (uo) is det.

%-----------------------------------------------------------------------------%

:- func (vector) / (vector) = (vector).
:- mode (in) / (in) = (uo) is det.
:- mode (di) / (in) = (uo) is det.
:- mode (in) / (di) = (uo) is det.
:- mode (di) / (di) = (uo) is det.

%-----------------------------------------------------------------------------%

:- func -(vector) = (vector).
:- mode -(di) = (uo) is det.
:- mode -(in) = (uo) is det.

%-----------------------------------------------------------------------------%

:- func scale(vector, float) = (vector) .
:- mode scale(in, in) = (uo) is det.
:- mode scale(di, in) = (uo) is det.

%-----------------------------------------------------------------------------%

:- func dot(vector::in, vector::in) = (float::uo) is det.

%=============================================================================%
:- implementation.
%=============================================================================%

:- import_module float.

%=============================================================================%
% Mercury implementation. This is used for all backends except C.
%=============================================================================%

:- type vector ---> vec(vx::float, vy::float, vz::float, vw::float).

%-----------------------------------------------------------------------------%

vector(X::in, Y::in, Z::in, W::in) = (vec(X+0.0, Y+0.0, Z+0.0, W+0.0)::uo).
vector(X+0.0::uo, Y+0.0::uo, Z+0.0::uo, W+0.0::uo) = (vec(X, Y, Z, W)::in).

% Needed because of the disjunct bodies, even in pure Mercury.
:- pragma promise_pure(vector/4).

%-----------------------------------------------------------------------------%

x(V) = V ^ vx.
y(V) = V ^ vy.
z(V) = V ^ vz.
w(V) = V ^ vw.

%-----------------------------------------------------------------------------%

'x :='(V, X) = V ^ vx := X.
'y :='(V, Y) = V ^ vy := Y.
'z :='(V, Z) = V ^ vz := Z.
'w :='(V, W) = V ^ vw := W.

%-----------------------------------------------------------------------------%

vector(X1, Y1, Z1, W1) + vector(X2, Y2, Z2, W2) = vector(X1+X2, Y1+Y2, Z1+Z2, W1+W2).

%-----------------------------------------------------------------------------%

vector(X1, Y1, Z1, W1) - vector(X2, Y2, Z2, W2) = vector(X1-X2, Y1-Y2, Z1-Z2, W1-W2).

%-----------------------------------------------------------------------------%

vector(X1, Y1, Z1, W1) * vector(X2, Y2, Z2, W2) = vector(X1*X2, Y1*Y2, Z1*Z2, W1*W2).

%-----------------------------------------------------------------------------%

vector(X1, Y1, Z1, W1) / vector(X2, Y2, Z2, W2) = vector(X1/X2, Y1/Y2, Z1/Z2, W1/W2).

%-----------------------------------------------------------------------------%

-vector(X, Y, Z, W) = vector(-X, -Y, -Z, -W).

%-----------------------------------------------------------------------------%

scale(vector(X, Y, Z, W), T) = vector(X*T, Y*T, Z*T, W*T).

%-----------------------------------------------------------------------------%

dot(vector(X1, Y1, Z1, W1), vector(X2, Y2, Z2, W2)) = (X1*X2) + (Y1*Y2) + (Z1*Z2) + (W1*W2).

%=============================================================================%
% DirectXMath implementation, used in C backends.
%=============================================================================%

:- pragma foreign_decl("C", "#include ""dxmath.h"" ").
:- pragma foreign_decl("C", " struct MerDX_VECTOR4_wrapper; ").
:- pragma foreign_type("C", vector, "struct MerDX_VECTOR4_wrapper*").

:- pragma foreign_decl("C", 
    "
#define MerDX_AllocateVector4() \\
    ((struct MerDX_VECTOR4_wrapper*)MR_GC_malloc_atomic(sizeof(DXMATH_VECTOR4) + 15))
    
typedef DXMATH_VECTOR4* MerDX_DXMATH_VECTOR4_PTR;
    MerDX_DXMATH_VECTOR4_PTR DXMATH_CONSTEXPR DXMATH_CALL
    MerDX_Vector4Offset(struct MerDX_VECTOR4_wrapper *w){
        MR_Word data = (MR_Word)w;
        if((data & 0xF) == 0) return (DXMATH_VECTOR4*)w;
        data = (data-1) & (~((MR_Word)0xF));
        data += 0x10;
        return (DXMATH_VECTOR4*)data;
    }

typedef struct MerDX_VECTOR4_wrapper* MerDX_VECTOR4_PTR;
    MerDX_VECTOR4_PTR DXMATH_CONSTEXPR DXMATH_CALL
    MerDX_CreateVector4Assign(DXMATH_VECTOR4 V){
        struct MerDX_VECTOR4_wrapper *const wrapper = MerDX_AllocateVector4();
        *MerDX_Vector4Offset(wrapper) = V;
        return wrapper;
    }
    
#define MerDX_CreateVector4(FX, FY, FZ, FW) \\
    MerDX_CreateVector4Assign(DXMATH_Vector4Set(FX, FY, FZ, FW))
    
#define MerDX_Vector4Unwrap(VEC_WRAP) \\
    (*MerDX_Vector4Offset((VEC_WRAP)))

#define MerDX_DESTRUCTIVE_BINOP(MERDX_OP, MERDX_OVERWRITE, MERDX_V1, MERDX_V2) \\
    (MERDX_OVERWRITE); \\
    do{ \\
    const DXMATH_VECTOR4 MerDX_DESTRUCTIVE_BINOP_vec = \\
        MERDX_OP(MerDX_Vector4Unwrap(MERDX_V1), MerDX_Vector4Unwrap(MERDX_V2)); \\
    MerDX_Vector4Unwrap(Out) = MerDX_DESTRUCTIVE_BINOP_vec; \\
    } while(0)

    ").

%-----------------------------------------------------------------------------%

:- pragma foreign_proc("C", x(V::in) = (F::out),
    [will_not_call_mercury, will_not_throw_exception, promise_pure, thread_safe,
     will_not_modify_trail, does_not_affect_liveness, may_duplicate],
    " F = DXMATH_VectorChannel(MerDX_Vector4Unwrap(V), 0); ").

:- pragma foreign_proc("C", y(V::in) = (F::out),
    [will_not_call_mercury, will_not_throw_exception, promise_pure, thread_safe,
     will_not_modify_trail, does_not_affect_liveness, may_duplicate],
    " F = DXMATH_VectorChannel(MerDX_Vector4Unwrap(V), 1); ").

:- pragma foreign_proc("C", z(V::in) = (F::out),
    [will_not_call_mercury, will_not_throw_exception, promise_pure, thread_safe,
     will_not_modify_trail, does_not_affect_liveness, may_duplicate],
    " F = DXMATH_VectorChannel(MerDX_Vector4Unwrap(V), 2); ").

:- pragma foreign_proc("C", w(V::in) = (F::out),
    [will_not_call_mercury, will_not_throw_exception, promise_pure, thread_safe,
     will_not_modify_trail, does_not_affect_liveness, may_duplicate],
    " F = DXMATH_VectorChannel(MerDX_Vector4Unwrap(V), 3); ").
    
%-----------------------------------------------------------------------------%

:- pragma foreign_proc("C", 'x :='(V::in, F::in) = (Out::out),
    [will_not_call_mercury, will_not_throw_exception, promise_pure, thread_safe,
     will_not_modify_trail, does_not_affect_liveness, may_duplicate],
    "
    const DXMATH_VECTOR4 vec = MerDX_Vector4Unwrap(V);
    Out = MerDX_CreateVector4(
        F,
        DXMATH_VectorChannel(vec, 1),
        DXMATH_VectorChannel(vec, 2),
        DXMATH_VectorChannel(vec, 3)); ").

:- pragma foreign_proc("C", 'y :='(V::in, F::in) = (Out::out),
    [will_not_call_mercury, will_not_throw_exception, promise_pure, thread_safe,
     will_not_modify_trail, does_not_affect_liveness, may_duplicate],
    "
    const DXMATH_VECTOR4 vec = MerDX_Vector4Unwrap(V);
    Out = MerDX_CreateVector4(
        DXMATH_VectorChannel(vec, 0),
        F,
        DXMATH_VectorChannel(vec, 2),
        DXMATH_VectorChannel(vec, 3)); ").

:- pragma foreign_proc("C", 'z :='(V::in, F::in) = (Out::out),
    [will_not_call_mercury, will_not_throw_exception, promise_pure, thread_safe,
     will_not_modify_trail, does_not_affect_liveness, may_duplicate],
    "
    const DXMATH_VECTOR4 vec = MerDX_Vector4Unwrap(V);
    Out = MerDX_CreateVector4(
        DXMATH_VectorChannel(vec, 0),
        DXMATH_VectorChannel(vec, 1),
        F,
        DXMATH_VectorChannel(vec, 3)); ").

:- pragma foreign_proc("C", 'w :='(V::in, F::in) = (Out::out),
    [will_not_call_mercury, will_not_throw_exception, promise_pure, thread_safe,
     will_not_modify_trail, does_not_affect_liveness, may_duplicate],
    "
    const DXMATH_VECTOR4 vec = MerDX_Vector4Unwrap(V);
    Out = MerDX_CreateVector4(
        DXMATH_VectorChannel(vec, 0),
        DXMATH_VectorChannel(vec, 1),
        DXMATH_VectorChannel(vec, 2),
        F); ").

%-----------------------------------------------------------------------------%

:- pragma foreign_proc("C", vector(X::in, Y::in, Z::in, W::in) = (V::uo),
    [will_not_call_mercury, will_not_throw_exception, promise_pure, thread_safe,
     will_not_modify_trail, does_not_affect_liveness, may_duplicate],
    " V = MerDX_CreateVector4(X, Y, Z, W); ").

:- pragma foreign_proc("C", vector(X::uo, Y::uo, Z::uo, W::uo) = (V::in),
    [will_not_call_mercury, will_not_throw_exception, promise_pure, thread_safe,
     will_not_modify_trail, does_not_affect_liveness, may_duplicate],
    "
    const DXMATH_VECTOR4 vec = MerDX_Vector4Unwrap(V);
    X = DXMATH_VectorChannel(vec, 0);
    Y = DXMATH_VectorChannel(vec, 1);
    Z = DXMATH_VectorChannel(vec, 2);
    W = DXMATH_VectorChannel(vec, 3); ").

%-----------------------------------------------------------------------------%

:- pragma foreign_proc("C", (V1::in) + (V2::in) = (Out::uo),
    [will_not_call_mercury, will_not_throw_exception, promise_pure, thread_safe,
     will_not_modify_trail, does_not_affect_liveness, may_duplicate],
    "
    const DXMATH_VECTOR4 vec = DXMATH_Vector4Add(MerDX_Vector4Unwrap(V1), MerDX_Vector4Unwrap(V2));
    Out = MerDX_CreateVector4Assign(vec);
    ").

%-----------------------------------------------------------------------------%

:- pragma foreign_proc("C", (V1::di) + (V2::in) = (Out::uo),
    [will_not_call_mercury, will_not_throw_exception, promise_pure, thread_safe,
     will_not_modify_trail, does_not_affect_liveness, may_duplicate],
    " Out = MerDX_DESTRUCTIVE_BINOP(DXMATH_Vector4Add, V1, V1, V2); ").

:- pragma foreign_proc("C", (V1::in) + (V2::di) = (Out::uo),
    [will_not_call_mercury, will_not_throw_exception, promise_pure, thread_safe,
     will_not_modify_trail, does_not_affect_liveness, may_duplicate],
    " Out = MerDX_DESTRUCTIVE_BINOP(DXMATH_Vector4Add, V2, V1, V2); ").

:- pragma foreign_proc("C", (V1::di) + (V2::di) = (Out::uo),
    [will_not_call_mercury, will_not_throw_exception, promise_pure, thread_safe,
     will_not_modify_trail, does_not_affect_liveness, may_duplicate],
    " Out = MerDX_DESTRUCTIVE_BINOP(DXMATH_Vector4Add, V1, V1, V2); ").

:- pragma promise_pure(('+')/2).

%-----------------------------------------------------------------------------%

:- pragma foreign_proc("C", (V1::in) - (V2::in) = (Out::uo),
    [will_not_call_mercury, will_not_throw_exception, promise_pure, thread_safe,
     will_not_modify_trail, does_not_affect_liveness, may_duplicate],
    "
    const DXMATH_VECTOR4 vec = DXMATH_Vector4Subtract(MerDX_Vector4Unwrap(V1), MerDX_Vector4Unwrap(V2));
    Out = MerDX_CreateVector4Assign(vec);
    ").

:- pragma foreign_proc("C", (V1::di) - (V2::in) = (Out::uo),
    [will_not_call_mercury, will_not_throw_exception, promise_pure, thread_safe,
     will_not_modify_trail, does_not_affect_liveness, may_duplicate],
    " Out = MerDX_DESTRUCTIVE_BINOP(DXMATH_Vector4Subtract, V1, V1, V2); ").

:- pragma foreign_proc("C", (V1::in) - (V2::di) = (Out::uo),
    [will_not_call_mercury, will_not_throw_exception, promise_pure, thread_safe,
     will_not_modify_trail, does_not_affect_liveness, may_duplicate],
    " Out = MerDX_DESTRUCTIVE_BINOP(DXMATH_Vector4Subtract, V2, V1, V2); ").

:- pragma foreign_proc("C", (V1::di) - (V2::di) = (Out::uo),
    [will_not_call_mercury, will_not_throw_exception, promise_pure, thread_safe,
     will_not_modify_trail, does_not_affect_liveness, may_duplicate],
    " Out = MerDX_DESTRUCTIVE_BINOP(DXMATH_Vector4Subtract, V1, V1, V2); ").

:- pragma promise_pure(('-')/2).

%-----------------------------------------------------------------------------%

:- pragma foreign_proc("C", (V1::in) * (V2::in) = (Out::uo),
    [will_not_call_mercury, will_not_throw_exception, promise_pure, thread_safe,
     will_not_modify_trail, does_not_affect_liveness, may_duplicate],
    "
    const DXMATH_VECTOR4 vec = DXMATH_Vector4Multiply(MerDX_Vector4Unwrap(V1), MerDX_Vector4Unwrap(V2));
    Out = MerDX_CreateVector4Assign(vec);
    ").

:- pragma foreign_proc("C", (V1::di) * (V2::in) = (Out::uo),
    [will_not_call_mercury, will_not_throw_exception, promise_pure, thread_safe,
     will_not_modify_trail, does_not_affect_liveness, may_duplicate],
    " Out = MerDX_DESTRUCTIVE_BINOP(DXMATH_Vector4Multiply, V1, V1, V2); ").

:- pragma foreign_proc("C", (V1::in) * (V2::di) = (Out::uo),
    [will_not_call_mercury, will_not_throw_exception, promise_pure, thread_safe,
     will_not_modify_trail, does_not_affect_liveness, may_duplicate],
    " Out = MerDX_DESTRUCTIVE_BINOP(DXMATH_Vector4Multiply, V2, V1, V2); ").

:- pragma foreign_proc("C", (V1::di) * (V2::di) = (Out::uo),
    [will_not_call_mercury, will_not_throw_exception, promise_pure, thread_safe,
     will_not_modify_trail, does_not_affect_liveness, may_duplicate],
    " Out = MerDX_DESTRUCTIVE_BINOP(DXMATH_Vector4Multiply, V1, V1, V2); ").

:- pragma promise_pure(('*')/2).

%-----------------------------------------------------------------------------%

:- pragma foreign_proc("C", (V1::in) / (V2::in) = (Out::uo),
    [will_not_call_mercury, will_not_throw_exception, promise_pure, thread_safe,
     will_not_modify_trail, does_not_affect_liveness, may_duplicate],
    "
    const DXMATH_VECTOR4 vec = DXMATH_Vector4Divide(MerDX_Vector4Unwrap(V1), MerDX_Vector4Unwrap(V2));
    Out = MerDX_CreateVector4Assign(vec);
    ").

:- pragma foreign_proc("C", (V1::di) / (V2::in) = (Out::uo),
    [will_not_call_mercury, will_not_throw_exception, promise_pure, thread_safe,
     will_not_modify_trail, does_not_affect_liveness, may_duplicate],
    " Out = MerDX_DESTRUCTIVE_BINOP(DXMATH_Vector4Divide, V1, V1, V2); ").

:- pragma foreign_proc("C", (V1::in) / (V2::di) = (Out::uo),
    [will_not_call_mercury, will_not_throw_exception, promise_pure, thread_safe,
     will_not_modify_trail, does_not_affect_liveness, may_duplicate],
    " Out = MerDX_DESTRUCTIVE_BINOP(DXMATH_Vector4Divide, V2, V1, V2); ").

:- pragma foreign_proc("C", (V1::di) / (V2::di) = (Out::uo),
    [will_not_call_mercury, will_not_throw_exception, promise_pure, thread_safe,
     will_not_modify_trail, does_not_affect_liveness, may_duplicate],
    " Out = MerDX_DESTRUCTIVE_BINOP(DXMATH_Vector4Divide, V1, V1, V2); ").

:- pragma promise_pure(('/')/2).

%-----------------------------------------------------------------------------%

:- pragma foreign_proc("C", -(V::in) = (Out::uo),
    [will_not_call_mercury, will_not_throw_exception, promise_pure, thread_safe,
     will_not_modify_trail, does_not_affect_liveness, may_duplicate],
    "
    const DXMATH_VECTOR4 vec = DXMATH_Vector4Subtract(DXMATH_Vector4Zero(), MerDX_Vector4Unwrap(V));
    Out = MerDX_CreateVector4Assign(vec); ").

%-----------------------------------------------------------------------------%

:- pragma foreign_proc("C", -(V::di) = (Out::uo),
    [will_not_call_mercury, will_not_throw_exception, promise_pure, thread_safe,
     will_not_modify_trail, does_not_affect_liveness, may_duplicate],
    ";
    Out = V;
    const DXMATH_VECTOR4 vec = DXMATH_Vector4Subtract(DXMATH_Vector4Zero(), MerDX_Vector4Unwrap(V));
    MerDX_Vector4Unwrap(Out) = vec;
    ").
    
:- pragma promise_pure(('-')/1).

%-----------------------------------------------------------------------------%

:- pragma foreign_proc("C", scale(V::in, MerF::in) = (Out::uo),
    [will_not_call_mercury, will_not_throw_exception, promise_pure, thread_safe,
     will_not_modify_trail, does_not_affect_liveness, may_duplicate],
    "
    const float F = MerF;
    const DXMATH_VECTOR4 vec = DXMATH_Vector4Multiply(DXMATH_Vector4Set(F, F, F, F), MerDX_Vector4Unwrap(V));
    Out = MerDX_CreateVector4Assign(vec);
    ").

%-----------------------------------------------------------------------------%

:- pragma foreign_proc("C", scale(V::di, MerF::in) = (Out::uo),
    [will_not_call_mercury, will_not_throw_exception, promise_pure, thread_safe,
     will_not_modify_trail, does_not_affect_liveness, may_duplicate],
    "
    const float F = MerF;
    Out = V;
    const DXMATH_VECTOR4 vec = DXMATH_Vector4Multiply(DXMATH_Vector4Set(F, F, F, F), MerDX_Vector4Unwrap(V));
    MerDX_Vector4Unwrap(Out) = vec;
    ").

:- pragma promise_pure(scale/2).

%-----------------------------------------------------------------------------%

:- pragma foreign_proc("C", dot(V1::in, V2::in) = (Out::uo),
    [will_not_call_mercury, will_not_throw_exception, promise_pure, thread_safe,
     will_not_modify_trail, does_not_affect_liveness, may_duplicate],
    " 
    Out = DXMATH_Vector4Dot(MerDX_Vector4Unwrap(V1), MerDX_Vector4Unwrap(V2));
    ").
