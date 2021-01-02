/**
RADAU solver implementation, adapted from  
https://github.com/scipy/scipy/blob/master/scipy/integrate/_ivp/radau.py
*/

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


/**
Predict by which factor to increase/decrease the step size.
The algorithm is described in [1]_.
Parameters
----------
h_abs, h_abs_old : float
    Current and previous values of the step size, `h_abs_old` can be None
    (see Notes).
error_norm, error_norm_old : float
    Current and previous values of the error norm, `error_norm_old` can
    be None (see Notes).
Returns
-------
factor : float
    Predicted factor.
Notes
-----
If `h_abs_old` and `error_norm_old` are both not None then a two-step
algorithm is used, otherwise a one-step algorithm is used.
References
----------
.. [1] E. Hairer, S. P. Norsett G. Wanner, "Solving Ordinary Differential
        Equations II: Stiff and Differential-Algebraic Problems", Sec. IV.8.
*/
real predict_factor(real h_abs, real h_abs_old, real error_norm, real error_norm_old):
{
    if (error_norm_old is null || h_abs_old is null || error_norm == 0){
        auto multiplier = 1;
    }
    else{
        auto multiplier = h_abs / h_abs_old * (error_norm_old / error_norm) ^^ 0.25;
    }
    // with np.errstate(divide='ignore'):
    auto factor = min(1, multiplier) * error_norm ^^ -0.25;

    return factor;
}

/**
Implicit Runge-Kutta method of Radau IIA family of order 5.
The implementation follows [1]_. The error is controlled with a
third-order accurate embedded formula. A cubic polynomial which satisfies
the collocation conditions is used for the dense output.
Parameters
----------
fun : callable
    Right-hand side of the system. The calling signature is ``fun(t, y)``.
    Here ``t`` is a scalar, and there are two options for the ndarray ``y``:
    It can either have shape (n,); then ``fun`` must return array_like with
    shape (n,). Alternatively it can have shape (n, k); then ``fun``
    must return an array_like with shape (n, k), i.e., each column
    corresponds to a single column in ``y``. The choice between the two
    options is determined by `vectorized` argument (see below). The
    vectorized implementation allows a faster approximation of the Jacobian
    by finite differences (required for this solver).
t0 : float
    Initial time.
y0 : array_like, shape (n,)
    Initial state.
t_bound : float
    Boundary time - the integration won't continue beyond it. It also
    determines the direction of the integration.
first_step : float or None, optional
    Initial step size. Default is ``None`` which means that the algorithm
    should choose.
max_step : float, optional
    Maximum allowed step size. Default is np.inf, i.e., the step size is not
    bounded and determined solely by the solver.
rtol, atol : float and array_like, optional
    Relative and absolute tolerances. The solver keeps the local error
    estimates less than ``atol + rtol * abs(y)``. Here `rtol` controls a
    relative accuracy (number of correct digits). But if a component of `y`
    is approximately below `atol`, the error only needs to fall within
    the same `atol` threshold, and the number of correct digits is not
    guaranteed. If components of y have different scales, it might be
    beneficial to set different `atol` values for different components by
    passing array_like with shape (n,) for `atol`. Default values are
    1e-3 for `rtol` and 1e-6 for `atol`.
jac : {None, array_like, sparse_matrix, callable}, optional
    Jacobian matrix of the right-hand side of the system with respect to
    y, required by this method. The Jacobian matrix has shape (n, n) and
    its element (i, j) is equal to ``d f_i / d y_j``.
    There are three ways to define the Jacobian:
        * If array_like or sparse_matrix, the Jacobian is assumed to
            be constant.
        * If callable, the Jacobian is assumed to depend on both
            t and y; it will be called as ``jac(t, y)`` as necessary.
            For the 'Radau' and 'BDF' methods, the return value might be a
            sparse matrix.
        * If None (default), the Jacobian will be approximated by
            finite differences.
    It is generally recommended to provide the Jacobian rather than
    relying on a finite-difference approximation.
jac_sparsity : {None, array_like, sparse matrix}, optional
    Defines a sparsity structure of the Jacobian matrix for a
    finite-difference approximation. Its shape must be (n, n). This argument
    is ignored if `jac` is not `None`. If the Jacobian has only few non-zero
    elements in *each* row, providing the sparsity structure will greatly
    speed up the computations [2]_. A zero entry means that a corresponding
    element in the Jacobian is always zero. If None (default), the Jacobian
    is assumed to be dense.
vectorized : bool, optional
    Whether `fun` is implemented in a vectorized fashion. Default is False.
Attributes
----------
n : int
    Number of equations.
status : string
    Current status of the solver: 'running', 'finished' or 'failed'.
t_bound : float
    Boundary time.
direction : float
    Integration direction: +1 or -1.
t : float
    Current time.
y : ndarray
    Current state.
t_old : float
    Previous time. None if no steps were made yet.
step_size : float
    Size of the last successful step. None if no steps were made yet.
nfev : int
    Number of evaluations of the right-hand side.
njev : int
    Number of evaluations of the Jacobian.
nlu : int
    Number of LU decompositions.
References
----------
.. [1] E. Hairer, G. Wanner, "Solving Ordinary Differential Equations II:
        Stiff and Differential-Algebraic Problems", Sec. IV.8.
.. [2] A. Curtis, M. J. D. Powell, and J. Reid, "On the estimation of
        sparse Jacobian matrices", Journal of the Institute of Mathematics
        and its Applications, 13, pp. 117-120, 1974.
*/
class Radau{

}

//Work in progress

