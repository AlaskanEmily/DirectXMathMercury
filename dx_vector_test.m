% Any copyright is dedicated to the Public Domain.
% http://creativecommons.org/publicdomain/zero/1.0/

:- module dx_vector_test.
%=============================================================================%
:- interface.
%=============================================================================%

:- import_module dx_vector.
:- use_module io.

:- func frobulate(vector::di, vector::di, vector::di) = (vector::uo) is det.

:- pred main(io.io::di, io.io::uo) is det.

%=============================================================================%
:- implementation.
%=============================================================================%

frobulate(V1, V2, V3) =
    scale(vector(1.0, 1.0, 1.0, 1.0) / -V2, dot(V1 + V3, V1 - V3)).

main(!IO) :-
    Result = frobulate(
        vector(0.728361, 293.0, 0.001, 1.0),
        vector(777.0, 776.0, 778.0, 1.0),
        vector(0.9, 0.09, 216.1, 1.0)),
    
    io.write_string("X: ", !IO),
    io.write_float(Result ^ x, !IO), io.nl(!IO),
    io.write_string("Y: ", !IO),
    io.write_float(Result ^ y, !IO), io.nl(!IO),
    io.write_string("Z: ", !IO),
    io.write_float(Result ^ z, !IO), io.nl(!IO).
