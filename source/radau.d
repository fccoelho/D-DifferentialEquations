module radau;

import std.stdio;
import std.range;
import std.typecons : tuple, Tuple;
import mir.interpolate;
import mir.ndslice;
import std.typecons : tuple, Tuple;

const real S6 = 6 ^^ 0.5;

/// Butcher tableau. A is not used directly, see below.
real[] C = [(4 - S6) / 10, (4 + S6) / 10, 1];
real[] E = [(-13 - 7 * S6) / 3, (-13 + 7 * S6) / 3, -1 / 3];

/** Eigendecomposition of A is done: A = T L T^^-1. There is 1 real eigenvalue 
and a complex conjugate pair. They are written below.
 */
real MU_REAL = 3 + 3 ^^ (2 / 3) - 3 ^^ (1 / 3);
creal MU_COMPLEX = (3 + 0.5 * (3 ^^ (1 / 3) - 3 ^^ (2 / 3)) - 0.5i * (3 ^^ (5 / 6) + 3 ^^ (7 / 6)));

/// These are transformation matrices.
real[][] T = [
    [0.09443876248897524, -0.14125529502095421, 0.03002919410514742],
    [0.25021312296533332, 0.20412935229379994, -0.38294211275726192], [1, 1, 0]
];

enum real[][] TI = [
        [4.17871859155190428, 0.32768282076106237, 0.52337644549944951],
        [-4.17871859155190428, -0.32768282076106237, 0.47662355450055044],
        [0.50287263494578682, -2.57192694985560522, 0.59603920482822492]
    ];

/// These linear combinations are used in the algorithm.
auto TI_REAL = TI[0];
// creal TI_COMPLEX = TI[1] + 1i * TI[2];
creal[] TI_COMPLEX = [
    -4.17871859155190428 + 0.50287263494578682i,
    -0.32768282076106237 - 2.57192694985560522i,
    0.47662355450055044 + 0.59603920482822492i
];



/// Interpolator coefficients.
real[][] P = [
    [13 / 3 + 7 * S6 / 3, -23 / 3 - 22 * S6 / 3, 10 / 3 + 5 * S6],
    [13 / 3 - 7 * S6 / 3, -23 / 3 + 22 * S6 / 3, 10 / 3 - 5 * S6],
    [1 / 3, -8 / 3, 10 / 3]
];

uint NEWTON_MAXITER = 6; /// Maximum number of Newton iterations.
real MIN_FACTOR = 0.2; /// Minimum allowed decrease in a step size.
real MAX_FACTOR = 10; /// Maximum allowed increase in a step size.

//Work in progress
