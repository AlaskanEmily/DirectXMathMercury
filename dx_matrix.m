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

:- module dx_matrix.
%=============================================================================%
% Bindings to DirectX Math matrix 4x4.
% Provides pure Mercury fallbacks for grades that do not interface with C.
:- interface.
%=============================================================================%

:- use_module dx_vector.

:- type matrix. % ---> matrix(
%   ax::float, ay::float, az::float, aw::float,
%   bx::float, by::float, bz::float, bw::float,
%   cx::float, cy::float, cz::float, cw::float,
%   dx::float, dy::float, dz::float, dw::float).

%-----------------------------------------------------------------------------%

% Used to simulate a matrix functor.
:- func matrix(float, float, float, float,
    float, float, float, float,
    float, float, float, float,
    float, float, float, float) = matrix.
:- mode matrix(in, in, in, in,
    in, in, in, in,
    in, in, in, in,
    in, in, in, in) = (uo) is det.
:- mode matrix(out, out, out, out,
    out, out, out, out,
    out, out, out, out,
    out, out, out, out) = (in) is det.

%-----------------------------------------------------------------------------%

:- func matrix(dx_vector.vector, dx_vector.vector, dx_vector.vector, dx_vector.vector) = matrix.
:- mode matrix(in, in, in, in) = (uo) is det.
:- mode matrix(di, di, di, di) = (uo) is det.
:- mode matrix(uo, uo, uo, uo) = (in) is det.

%-----------------------------------------------------------------------------%

:- func a(matrix) = dx_vector.vector.
:- func b(matrix) = dx_vector.vector.
:- func c(matrix) = dx_vector.vector.
:- func d(matrix) = dx_vector.vector.

% Field access functions.
:- func ax(matrix) = float.
:- func ay(matrix) = float.
:- func az(matrix) = float.
:- func aw(matrix) = float.

:- func bx(matrix) = float.
:- func by(matrix) = float.
:- func bz(matrix) = float.
:- func bw(matrix) = float.

:- func cx(matrix) = float.
:- func cy(matrix) = float.
:- func cz(matrix) = float.
:- func cw(matrix) = float.

:- func dx(matrix) = float.
:- func dy(matrix) = float.
:- func dz(matrix) = float.
:- func dw(matrix) = float.

%-----------------------------------------------------------------------------%

:- func 'a :='(matrix, dx_vector.vector) = matrix.
:- func 'b :='(matrix, dx_vector.vector) = matrix.
:- func 'c :='(matrix, dx_vector.vector) = matrix.
:- func 'd :='(matrix, dx_vector.vector) = matrix.

:- func 'ax :='(matrix, float) = matrix.
:- func 'ay :='(matrix, float) = matrix.
:- func 'az :='(matrix, float) = matrix.
:- func 'aw :='(matrix, float) = matrix.

:- func 'bx :='(matrix, float) = matrix.
:- func 'by :='(matrix, float) = matrix.
:- func 'bz :='(matrix, float) = matrix.
:- func 'bw :='(matrix, float) = matrix.

:- func 'cx :='(matrix, float) = matrix.
:- func 'cy :='(matrix, float) = matrix.
:- func 'cz :='(matrix, float) = matrix.
:- func 'cw :='(matrix, float) = matrix.

:- func 'dx :='(matrix, float) = matrix.
:- func 'dy :='(matrix, float) = matrix.
:- func 'dz :='(matrix, float) = matrix.
:- func 'dw :='(matrix, float) = matrix.

:- func transpose(matrix) = matrix.
:- mode transpose(in) = (uo) is det.
:- mode transpose(di) = (uo) is det.

:- pred transpose(matrix, matrix).
:- mode transpose(in, uo) is det.
:- mode transpose(di, uo) is det.
:- mode transpose(uo, in) is det.
:- mode transpose(uo, di) is det.

:- func identity = matrix.

%=============================================================================%
:- implementation.
%=============================================================================%

% Never learn how the sausage is made.

:- import_module float.
:- use_module math.

% Avoid allocating a new matrix just for identity. This is a nop in pure Mercury.
:- impure pred initalize_identity is det.
:- initialise initalize_identity/0.

% We try to force inlining as much as possible. This allows many equations to
% actually compile down from Mercury all the way to ASM which is pretty close
% to what the C version would have been (sans some extra allocations ;_; )

:- pragma inline(a/1).
:- pragma inline(b/1).
:- pragma inline(c/1).
:- pragma inline(d/1).

:- pragma inline(ax/1).
:- pragma inline(ay/1).
:- pragma inline(az/1).
:- pragma inline(aw/1).

:- pragma inline(bx/1).
:- pragma inline(by/1).
:- pragma inline(bz/1).
:- pragma inline(bw/1).

:- pragma inline(cx/1).
:- pragma inline(cy/1).
:- pragma inline(cz/1).
:- pragma inline(cw/1).

:- pragma inline(dx/1).
:- pragma inline(dy/1).
:- pragma inline(dz/1).
:- pragma inline(dw/1).

:- pragma inline('a :='/2).
:- pragma inline('b :='/2).
:- pragma inline('c :='/2).
:- pragma inline('d :='/2).

:- pragma inline('ax :='/2).
:- pragma inline('ay :='/2).
:- pragma inline('az :='/2).
:- pragma inline('aw :='/2).

:- pragma inline('bx :='/2).
:- pragma inline('by :='/2).
:- pragma inline('bz :='/2).
:- pragma inline('bw :='/2).

:- pragma inline('cx :='/2).
:- pragma inline('cy :='/2).
:- pragma inline('cz :='/2).
:- pragma inline('cw :='/2).

:- pragma inline('dx :='/2).
:- pragma inline('dy :='/2).
:- pragma inline('dz :='/2).
:- pragma inline('dw :='/2).

:- pragma inline(matrix/16).
:- pragma inline(matrix/4).
:- pragma inline(transpose/1).
:- pragma inline(transpose/2).
:- pragma inline(identity/0).

%=============================================================================%
% Pred and func versions of operators come first, as these are shared between
% both the C backend and the Mercury backend.
% This also means that the C and Mercury backends only need to define the
% vector field-access funcs, since both matrix types are implemented in terms
% of four dx_vector.vector's.
%=============================================================================%

:- func copy_vec(dx_vector.vector) = (dx_vector.vector).
:- mode copy_vec(in) = (uo) is det.
:- mode copy_vec(uo) = (in) is det.
copy_vec(dx_vector.vector(X, Y, Z, W)) = dx_vector.vector(X, Y, Z, W).

% This is what the C code would basically do anyway, and since the Mercury
% implementation uses vectors under the hood, this also works there.
matrix(
    AX, AY, AZ, AW,
    BX, BY, BZ, BW,
    CX, CY, CZ, CW,
    DX, DY, DZ, DW) = matrix(A, B, C, D) :-
    
    A = dx_vector.vector(AX, AY, AZ, AW),
    B = dx_vector.vector(BX, BY, BZ, BW),
    C = dx_vector.vector(CX, CY, CZ, CW),
    D = dx_vector.vector(DX, DY, DZ, DW).

% Destructuring vectors is pretty much needed here regardless.
ax(M) = M ^ a ^ dx_vector.x.
ay(M) = M ^ a ^ dx_vector.y.
az(M) = M ^ a ^ dx_vector.z.
aw(M) = M ^ a ^ dx_vector.w.

bx(M) = M ^ b ^ dx_vector.x.
by(M) = M ^ b ^ dx_vector.y.
bz(M) = M ^ b ^ dx_vector.z.
bw(M) = M ^ b ^ dx_vector.w.

cx(M) = M ^ c ^ dx_vector.x.
cy(M) = M ^ c ^ dx_vector.y.
cz(M) = M ^ c ^ dx_vector.z.
cw(M) = M ^ c ^ dx_vector.w.

dx(M) = M ^ d ^ dx_vector.x.
dy(M) = M ^ d ^ dx_vector.y.
dz(M) = M ^ d ^ dx_vector.z.
dw(M) = M ^ d ^ dx_vector.w.

% TODO: di/uo variants might be smart :)
'ax :='(M, X) = M ^ a := ((M ^ a) ^ dx_vector.x := X).
'ay :='(M, Y) = M ^ a := ((M ^ a) ^ dx_vector.y := Y).
'az :='(M, Z) = M ^ a := ((M ^ a) ^ dx_vector.z := Z).
'aw :='(M, W) = M ^ a := ((M ^ a) ^ dx_vector.w := W).

'bx :='(M, X) = M ^ b := ((M ^ b) ^ dx_vector.x := X).
'by :='(M, Y) = M ^ b := ((M ^ b) ^ dx_vector.y := Y).
'bz :='(M, Z) = M ^ b := ((M ^ b) ^ dx_vector.z := Z).
'bw :='(M, W) = M ^ b := ((M ^ b) ^ dx_vector.w := W).

'cx :='(M, X) = M ^ c := ((M ^ c) ^ dx_vector.x := X).
'cy :='(M, Y) = M ^ c := ((M ^ c) ^ dx_vector.y := Y).
'cz :='(M, Z) = M ^ c := ((M ^ c) ^ dx_vector.z := Z).
'cw :='(M, W) = M ^ c := ((M ^ c) ^ dx_vector.w := W).

'dx :='(M, X) = M ^ d := ((M ^ d) ^ dx_vector.x := X).
'dy :='(M, Y) = M ^ d := ((M ^ d) ^ dx_vector.y := Y).
'dz :='(M, Z) = M ^ d := ((M ^ d) ^ dx_vector.z := Z).
'dw :='(M, W) = M ^ d := ((M ^ d) ^ dx_vector.w := W).

transpose(A::in, B::uo) :- B = transpose(A).
transpose(A::di, B::uo) :- B = transpose(A).
transpose(A::uo, B::in) :- A = transpose(B).
transpose(A::uo, B::di) :- A = transpose(B).

:- pragma promise_pure(transpose/2).

%=============================================================================%
% Mercury implementation. This is used for all backends except C.
%=============================================================================%

initalize_identity.

:- type matrix ---> mat(
    ma::dx_vector.vector,
    mb::dx_vector.vector, 
    mc::dx_vector.vector, 
    md::dx_vector.vector).

matrix(A, B, C, D) = mat(copy_vec(A), copy_vec(B), copy_vec(C), copy_vec(D)).

a(M) = M ^ ma.
b(M) = M ^ mb.
c(M) = M ^ mc.
d(M) = M ^ md.

'a :='(M, V) = M ^ ma := V.
'b :='(M, V) = M ^ mb := V.
'c :='(M, V) = M ^ mc := V.
'd :='(M, V) = M ^ md := V.

identity = matrix(
    1.0, 0.0, 0.0, 0.0,
    0.0, 1.0, 0.0, 0.0,
    0.0, 0.0, 1.0, 0.0,
    0.0, 0.0, 0.0, 1.0).

%=============================================================================%
% DirectXMath implementation, used in C backends.
%=============================================================================%

:- pragma foreign_decl("C", " #include ""dx_math_mercury.h"" ").
:- pragma foreign_decl("C", " extern DXMATH_MATRIX44 MerDX_IdentityMatrix_; ").
:- pragma foreign_decl("C", " extern MerDX_MATRIX44_PTR MerDX_IdentityMatrix; ").

% Avoid allocating a new matrix just for identity.
:- pragma foreign_code("C", "
    DXMATH_MATRIX44 MerDX_IdentityMatrix_;
    MerDX_MATRIX44_PTR MerDX_IdentityMatrix = NULL;
").

:- pragma foreign_type("C", matrix, "struct MerDX_MATRIX44_wrapper*").

:- pragma foreign_export("C", a(in) = (out), "MerDX_Matrix44CopyA").
:- pragma foreign_export("C", b(in) = (out), "MerDX_Matrix44CopyB").
:- pragma foreign_export("C", c(in) = (out), "MerDX_Matrix44CopyC").
:- pragma foreign_export("C", d(in) = (out), "MerDX_Matrix44CopyD").

:- pragma foreign_proc("C", initalize_identity,
    [will_not_call_mercury, will_not_throw_exception],
    "

#ifdef DXMATH_NO_INTRINSICS
    const MR_Word OK = 1;
#else
    const MR_Word OK = (((MR_Word)(&MerDX_IdentityMatrix_)) & 0xF) == 0;
#endif

    if(MerDX_IdentityMatrix == NULL){
        if(OK){
            DXMATH_Matrix44IdentityOut(&MerDX_IdentityMatrix_);
            MerDX_IdentityMatrix = (MerDX_MATRIX44_PTR)&MerDX_IdentityMatrix_;
        }
        else{
            MerDX_MATRIX44_PTR mat = MerDX_AllocateMatrix44();
            MerDX_Matrix44Unwrap(mat) = DXMATH_Matrix44Identity();
            MerDX_IdentityMatrix = mat;
        }
    }
    ").

:- pragma foreign_proc("C", matrix(A::in, B::in, C::in, D::in) = (Out::uo),
    [will_not_call_mercury, will_not_throw_exception, promise_pure, thread_safe,
     will_not_modify_trail, does_not_affect_liveness, may_duplicate],
    "
    DXMATH_MATRIX44 mat;
    DXMATH_Matrix44Row(mat, 0) = MerDX_Vector4Unwrap(A);
    DXMATH_Matrix44Row(mat, 1) = MerDX_Vector4Unwrap(B);
    DXMATH_Matrix44Row(mat, 2) = MerDX_Vector4Unwrap(C);
    DXMATH_Matrix44Row(mat, 3) = MerDX_Vector4Unwrap(D);
    Out = MerDX_CreateMatrix44Assign(mat);
    ").

:- pragma foreign_proc("C", matrix(A::di, B::di, C::di, D::di) = (Out::uo),
    [will_not_call_mercury, will_not_throw_exception, promise_pure, thread_safe,
     will_not_modify_trail, does_not_affect_liveness, may_duplicate],
    "
    DXMATH_MATRIX44 mat;
    
    DXMATH_Matrix44Row(mat, 0) = MerDX_Vector4Unwrap(A);
    DXMATH_Matrix44Row(mat, 1) = MerDX_Vector4Unwrap(B);
    DXMATH_Matrix44Row(mat, 2) = MerDX_Vector4Unwrap(C);
    DXMATH_Matrix44Row(mat, 3) = MerDX_Vector4Unwrap(D);
    
    MR_GC_Free(A);
    MR_GC_Free(B);
    MR_GC_Free(C);
    MR_GC_Free(D);
    
    Out = MerDX_CreateMatrix44Assign(mat);
    ").

:- pragma foreign_proc("C", matrix(A::uo, B::uo, C::uo, D::uo) = (M::in),
    [will_not_call_mercury, will_not_throw_exception, promise_pure, thread_safe,
     will_not_modify_trail, does_not_affect_liveness, may_duplicate],
    "
    A = MerDX_Matrix44CopyA(M);
    B = MerDX_Matrix44CopyB(M);
    C = MerDX_Matrix44CopyC(M);
    D = MerDX_Matrix44CopyD(M);
    ").

:- pragma promise_pure(matrix/4).

%-----------------------------------------------------------------------------%

:- pragma foreign_proc("C", a(M::in) = (V::out),
    [will_not_call_mercury, will_not_throw_exception, promise_pure, thread_safe,
     will_not_modify_trail, does_not_affect_liveness, may_duplicate],
    "
    const DXMATH_VECTOR4 vec = DXMATH_Matrix44Row(MerDX_Matrix44Unwrap(M), 0);
    V = MerDX_CreateVector4Assign(vec);
    ").

:- pragma foreign_proc("C", b(M::in) = (V::out),
    [will_not_call_mercury, will_not_throw_exception, promise_pure, thread_safe,
     will_not_modify_trail, does_not_affect_liveness, may_duplicate],
    "
    const DXMATH_VECTOR4 vec = DXMATH_Matrix44Row(MerDX_Matrix44Unwrap(M), 1);
    V = MerDX_CreateVector4Assign(vec);
    ").

:- pragma foreign_proc("C", c(M::in) = (V::out),
    [will_not_call_mercury, will_not_throw_exception, promise_pure, thread_safe,
     will_not_modify_trail, does_not_affect_liveness, may_duplicate],
    "
    const DXMATH_VECTOR4 vec = DXMATH_Matrix44Row(MerDX_Matrix44Unwrap(M), 2);
    V = MerDX_CreateVector4Assign(vec);
    ").

:- pragma foreign_proc("C", d(M::in) = (V::out),
    [will_not_call_mercury, will_not_throw_exception, promise_pure, thread_safe,
     will_not_modify_trail, does_not_affect_liveness, may_duplicate],
    "
    const DXMATH_VECTOR4 vec = DXMATH_Matrix44Row(MerDX_Matrix44Unwrap(M), 3);
    V = MerDX_CreateVector4Assign(vec);
    ").

:- pragma foreign_proc("C", 'a :='(M::in, V::in) = (Out::out),
    [will_not_call_mercury, will_not_throw_exception, promise_pure, thread_safe,
     will_not_modify_trail, does_not_affect_liveness, may_duplicate],
    " 
    DXMATH_MATRIX44 mat = MerDX_Matrix44Unwrap(M);
    DXMATH_Matrix44Row(mat, 0) = MerDX_Vector4Unwrap(V);
    Out = MerDX_CreateMatrix44Assign(mat);
    ").

:- pragma foreign_proc("C", 'b :='(M::in, V::in) = (Out::out),
    [will_not_call_mercury, will_not_throw_exception, promise_pure, thread_safe,
     will_not_modify_trail, does_not_affect_liveness, may_duplicate],
    " 
    DXMATH_MATRIX44 mat = MerDX_Matrix44Unwrap(M);
    DXMATH_Matrix44Row(mat, 1) = MerDX_Vector4Unwrap(V);
    Out = MerDX_CreateMatrix44Assign(mat);
    ").

:- pragma foreign_proc("C", 'c :='(M::in, V::in) = (Out::out),
    [will_not_call_mercury, will_not_throw_exception, promise_pure, thread_safe,
     will_not_modify_trail, does_not_affect_liveness, may_duplicate],
    " 
    DXMATH_MATRIX44 mat = MerDX_Matrix44Unwrap(M);
    DXMATH_Matrix44Row(mat, 2) = MerDX_Vector4Unwrap(V);
    Out = MerDX_CreateMatrix44Assign(mat);
    ").

:- pragma foreign_proc("C", 'd :='(M::in, V::in) = (Out::out),
    [will_not_call_mercury, will_not_throw_exception, promise_pure, thread_safe,
     will_not_modify_trail, does_not_affect_liveness, may_duplicate],
    " 
    DXMATH_MATRIX44 mat = MerDX_Matrix44Unwrap(M);
    DXMATH_Matrix44Row(mat, 3) = MerDX_Vector4Unwrap(V);
    Out = MerDX_CreateMatrix44Assign(mat);
    ").

:- pragma foreign_proc("C", transpose(M::in) = (Out::uo),
    [will_not_call_mercury, will_not_throw_exception, promise_pure, thread_safe,
     will_not_modify_trail, does_not_affect_liveness, may_duplicate],
    " 
    const DXMATH_MATRIX44 mat = DXMATH_TransposeMatrix(MerDX_Matrix44Unwrap(M));
    Out = MerDX_CreateMatrix44Assign(mat);
    ").

:- pragma foreign_proc("C", transpose(M::di) = (Out::uo),
    [will_not_call_mercury, will_not_throw_exception, promise_pure, thread_safe,
     will_not_modify_trail, does_not_affect_liveness, may_duplicate],
    "
    Out = M;
    const DXMATH_MATRIX44 mat = DXMATH_TransposeMatrix(MerDX_Matrix44Unwrap(M));
    MerDX_Matrix44Unwrap(Out) = mat;
    ").

:- pragma foreign_proc("C", identity = (Out::out),
    [will_not_call_mercury, will_not_throw_exception, promise_pure, thread_safe,
     will_not_modify_trail, does_not_affect_liveness, may_duplicate],
    "
        Out = (MerDX_MATRIX44_PTR)MerDX_IdentityMatrix;
    ").

:- pragma promise_pure(transpose/1).

