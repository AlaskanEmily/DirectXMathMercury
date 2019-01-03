# Any copyright is dedicated to the Public Domain.
# http://creativecommons.org/publicdomain/zero/1.0/

all:
	mmc --make dx_vector_test --c-include-directory -O6 dxmath_c --cflag -DDXMATH_NO_INTRINSICS
