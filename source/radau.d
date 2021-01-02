module radau;

import std.stdio;
import std.range;
import std.complex;
import std.typecons : tuple, Tuple;
import mir.interpolate;
import mir.ndslice;
import kaleidic.lubeck: mtimes;
import std.typecons : tuple, Tuple;

const real S6 = 6 ^^ 0.5;

/// Butcher tableau. A is not used directly, see below.
real[] C = [(4 - S6) / 10, (4 + S6) / 10, 1];
real[] E = [(-13 - 7 * S6) / 3, (-13 + 7 * S6) / 3, -1 / 3];

/** Eigendecomposition of A is done: A = T L T^^-1. There is 1 real eigenvalue 
and a complex conjugate pair. They are written below.
 */
real MU_REAL = 3 + 3 ^^ (2 / 3) - 3 ^^ (1 / 3);
auto MU_COMPLEX = (3 + 0.5 * (3 ^^ (1 / 3) - 3 ^^ (2 / 3)) - 0.5i * (3 ^^ (5 / 6) + 3 ^^ (7 / 6)));

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
cdouble[] TI_COMPLEX = [
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

/// Model function signature: F(t, y, args)
alias mfun = real[]delegate(real, real[], real[]) pure nothrow @safe;

/**
Solve the collocation system.
    Parameters
    ----------
    fun : callable
        Right-hand side of the system.
    t : float
        Current time.
    y : ndarray, shape (n,)
        Current state.
    h : float
        Step to try.
    Z0 : ndarray, shape (3, n)
        Initial guess for the solution. It determines new values of `y` at
        ``t + h * C`` as ``y + Z0``, where ``C`` is the Radau method constants.
    scale : float
        Problem tolerance scale, i.e. ``rtol * abs(y) + atol``.
    tol : float
        Tolerance to which solve the system. This value is compared with
        the normalized by `scale` error.
    LU_real, LU_complex
        LU decompositions of the system Jacobians.
    solve_lu : callable
        Callable which solves a linear system given a LU decomposition. The
        signature is ``solve_lu(LU, b)``.
    Returns
    -------
    converged : bool
        Whether iterations converged.
    n_iter : int
        Number of completed iterations.
    Z : ndarray, shape (3, n)
        Found solution.
    rate : float
        The rate of convergence.
*/
Tuple!(bool,int, real[], real) solve_collocation_system(mfun fun, real t, real[] y, real h, real[][] Z0, real scale, real tol,
                             real[][] LU_real, creal[][] LU_complex)
{
    auto n = y.length;
    auto M_real = MU_REAL / h;
    auto M_complex = MU_COMPLEX / h;

    auto W = TI.mtimes(Z0);
    auto Z = Z0.dup;

    auto F = repeat(0, 3*n).array.sliced(3,n); //np.empty((3, n));
    auto ch = h * C;

    auto dW_norm_old = null;
    auto dW = W.dup; dW[] *= 0;
    bool converged = False;
    auto rate = null;
    foreach(k; NEWTON_MAXITER.iota){
        foreach(i; 3.iota){
            F[i] = fun(t + ch[i], y + Z[i]);
        }
        if (!all(F.map!isfinite))
            break;

        auto f_real = F.T.mtimes(TI_REAL) - M_real * W[0];
        auto f_complex = F.T.mtimes(TI_COMPLEX) - M_complex * (W[1] + 1i * W[2]);

        auto dW_real = luSolve(LU_real, f_real);
        auto dW_complex = luSolve(LU_complex, f_complex);

        dW[0] = dW_real;
        dW[1] = dW_complex.re;
        dW[2] = dW_complex.im;

        dW_norm = norm(dW / scale);
        if (dW_norm_old !is null){
            rate = dW_norm / dW_norm_old;
        }
        if (rate !is null && (rate >= 1 || rate ^^ (NEWTON_MAXITER - k) / (1 - rate) * dW_norm > tol))
            break;

        W += dW;
        Z = T.mtimes(W);

        if (dW_norm == 0 || rate !is null && rate / (1 - rate) * dW_norm < tol){
            converged = true;
            break;
                }

        dW_norm_old = dW_norm;
    }

    return tuple(converged, k + 1, Z, rate);
}

//Work in progress

