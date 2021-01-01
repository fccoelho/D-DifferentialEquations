import std.stdio;
import std.math;
import std.range;
import std.array : array;
import std.algorithm;
import std.exception;
import std.typecons : tuple, Tuple;
import std.datetime.stopwatch;
import std.string;

/// function signature: F(t, y)
alias mfun = real[]delegate(real, real[], real[]) pure nothrow @safe;

Tuple!(real[], real[][]) odeint(mfun f, real[] y0, real[] times, real[] args = [
        ], string method="rk4", real tol = 1e-6, real hmax = 0,
        real hmin = 1e-6, uint maxiter = 1000)
{
    Tuple!(real[], real[][]) res;
    if (hmax == 0){
        hmax = hmin;
    }

    if (method == "rk4")
    {
        res = rk4(f, y0, times, tol, args = args,);
    }
    else if (method == "dopri")
    {
        res = dopri(f, times[0], y0, times[$ - 1], tol = tol, hmax = hmax,
                hmin = hmin, maxiter = maxiter, args = args);
    }
    return res;
}

/**
Runge-kutta 4th order step
*/
real[] rk4Step(mfun f, real x, real[] y, real h, real[] k1, real[] args) pure nothrow @safe
{
    real[] k2, k3, k4;
    auto s1 = y.dup;
    s1[] += 0.5 * h * k1[]; //new state, after dt=h/2 
    k2 = f(x + 0.5 * h, s1, args); // new dy after taking the step
    auto s2 = y.dup;
    s2[] += 0.5 * h * k2[];
    k3 = f(x + 0.5 * h, s2, args);
    auto s3 = y.dup;
    s3[] += h * k3[];
    k4 = f(x + h, s3, args);
    auto dy = k1.dup;
    dy[] += (2 * (k2[] + k3[]) + k4[]);
    dy[] /= 6.0;
    return dy;

}
/**
Two step average.
*/
real[] rk4TwoStep(mfun f, real x, real[] y, real h, real[] k1, real[] args) pure nothrow @safe
{
    real x1;
    real[] y1, step1, step2;
    step1 = rk4Step(f, x, y, 0.5 * h, k1, args);
    x1 = x + 0.5 * h;
    y1 = y.dup;
    y1[] += 0.5 * h * step1[];
    step2 = rk4Step(f, x1, y1, 0.5 * h, f(x1, y1, args), args);
    auto mean = step1.dup;
    mean[] += step2[];
    mean[] /= 2.0;
    return mean;
}

/**
Runge-Kutta, 4th order solver
f - function f(t,Y,args) returning the derivatives at time t
y0 - initial values
times - times for the solutions
tol - error tolerance
args - parameters to be passed to f
*/
Tuple!(real[], real[][]) rk4(mfun f, real[] y0, real[] times, const real tol = 1e-6, real[] args = [
        ]) @safe
{
    real[][] yout;
    real[] tout;
    yout ~= y0;

    auto h = times[1] - times[0];
    auto immutable hmax = abs(times[$ - 1] - times[0]);
    real hmin = tol;

    auto ycurr = y0;
    auto ylast = y0;
    auto qcurr = y0;
    auto qlast = y0;
    auto tcurr = times[0];
    auto tlast = times[0];
    auto fcurr = f(tcurr, ycurr, args);
    auto flast = f(tcurr, ycurr, args);
    real totalerr, steperr, scale;
    real totalvar = 0.0;

    foreach (immutable i, t; times[1 .. $])
    {
        steperr = 2 * tol;
        while ((t - tcurr) * h > 0)
        {
            // writeln(steperr > tol, ' ', steperr-tol,' ',t);
            auto const k1 = rk4Step(f, tcurr, ycurr, h, fcurr, args);
            auto k2 = rk4TwoStep(f, tcurr, ycurr, h, fcurr, args);

            scale = k2.map!abs.maxElement;
            /// steperr = max(abs(k1-k2))/2
            auto err = k1.dup;
            err[] -= k2[]; // error = k1-k2 
            steperr = err.map!abs.maxElement / 2;
            // if (steperr>0){
            //     writefln("t: %s, Err: %s, tol: %s, tcurr: %s, h: %s",t,steperr, tol, tcurr, h);
            // }
            /// compute the ideal step size factor and sanitize the result to prevent ridiculous changes
            auto hfac = ((tol * scale) / (steperr)) ^^ 0.25;
            hfac = min(10, max(0.01, hfac));
            // writeln(t, ' ', hfac, " ", steperr);
            /// repeat the step if there is a significant step size correction
            if (abs(h * hfac) < hmax && ((0.6 > hfac) || (hfac > 3)))
            {
                // recompute with new step size
                h = max(hmin, h * hfac);
                k2 = rk4TwoStep(f, tcurr, ycurr, h, fcurr, args);
            }
            // update and cycle the integration points
            ylast = ycurr.dup;
            ycurr[] += h * k2[];
            if (ycurr.any!isNaN)
            {
                writeln(h, k2);
            }
            tlast = tcurr;
            tcurr += h;
            // writefln("ycurr: %s, tcurr: %s, h: %s", ycurr,tcurr,h);
            flast = fcurr.dup;
            fcurr = f(tcurr, ycurr, args);
            // cubic Bezier control points
            qlast = ylast.dup;
            qlast[] += (tcurr - tlast) / 3 * flast[];
            qcurr = ycurr.dup;
            qcurr[] -= (tcurr - tlast) / 3 * fcurr[];

            totalvar += h * scale;
            totalerr = (1 + h * scale) * totalerr + h * steperr;
            //auto reportstr = "internal step to t=%12.8f \t" % tcurr;
        }
        //now tlast <= t <= tcurr, can interpolate the value for yout[i+1] using the cubic Bezier formula
        // writefln("Accepted step: t: %s, y: %s, Err: %s", tcurr, ycurr, steperr);
        auto s = (t - tlast) / (tcurr - tlast);
        auto res = ylast.dup;
        res[] *= 0;
        res[] += (1 - s) ^^ 2 * ((1 - s) * ylast[] + 3 * s * qlast[]) + s ^^ 2 * (
                3 * (1 - s) * qcurr[] + s * ycurr[]);
        // if (res.any!isNaN){
        //         writefln("h: %s, s: %s, Err: %s, (tcurr-tlast): %s, (t-tlast): %s",h,s, steperr, tcurr-tlast,t-tlast);
        //         }
        yout ~= res;
        tout ~= tcurr;
    }
    return tuple(tout, yout);
}

unittest
{
    real epsilon = 0.1;
    real b = 2.0;
    real c = 0.2;
    // model function
    real[] fun(real t, real[] Y, real[] P)
    {
        auto u = Y[0];
        auto v = Y[1];
        real[] dy = [0, 0];
        dy[0] += u - u ^^ 3 - v;
        dy[1] += P[0] * (u - P[1] * v - P[2]);
        return dy;
    }

    real[] inits = [0.1, 0.2];
    real[] times;
    times ~= 0;
    foreach (i; 0 .. 1000)
    {
        times ~= times[$ - 1] + 0.1;
    }

    auto res = rk4(&fun, inits, times, 1e-6, [epsilon, b, c]);
    // writeln(res);

}

/**
Synopsis
	  real dopri(real[] fxy(double x, double y),
					  double x0, double y0, double x1, double tol,
					  double hmax,  double hmin, int maxiter)
 
Parameters
	  fxy			   Input: derivative function y' = f(x, y)
						   y is the dependent variable, x is the independent
						   variable
	  x0, y0			Input: initial points, x0 <= x <= x1	y(x0) = y0
	  x1				Input: final value of x
	  tol			   Input: tolerance
	  hmax			  Input: maximum step size
	  hmin			  Input: minimum step size
	  maxiter		   Input: maximum number of iterations
	  flag			  Input: return flag
						   0   no errors
						   1   hmin exceeded
						   2   maximum iterations exceeded
 
Return value
	  value of y at last step x
 
Description
	  The routine dopri() implements the Dormand-Prince method of
	  solving an ordinary differential equation of the first order
	  y' = f(x,y).
 
Reference
	  The coefficients were obtained from
 
		  E.Hairer, S.P.Norsett and G.Wanner[1991],
			 "Solving Differential Equations I, Nonstiff Problems",
			 2e, Springer-Verlag, p. 178
 
WARNING
	  Check the flag after calling this routine!
 
Revisions
	  1998.05.02	  first version
*/
Tuple!(real[], real[][]) dopri(mfun fxy, real x0, real[] y0, real x1, real tol,
        real hmax, real hmin, uint maxiter, real[] args = [])
{
    enum real a21 = (1.0 / 5.0);
    enum real a31 = (3.0 / 40.0);
    enum real a32 = (9.0 / 40.0);
    enum real a41 = (44.0 / 45.0);
    enum real a42 = (-56.0 / 15.0);
    enum real a43 = (32.0 / 9.0);
    enum real a51 = (19_372.0 / 6561.0);
    enum real a52 = (-25_360.0 / 2187.0);
    enum real a53 = (64448.0 / 6561.0);
    enum real a54 = (-212.0 / 729.0);
    enum real a61 = (9017.0 / 3168.0);
    enum real a62 = (-355.0 / 33.0);
    enum real a63 = (46732.0 / 5247.0);
    enum real a64 = (49.0 / 176.0);
    enum real a65 = (-5103.0 / 18656.0);
    enum real a71 = (35.0 / 384.0);
    enum real a72 = (0.0);
    enum real a73 = (500.0 / 1113.0);
    enum real a74 = (125.0 / 192.0);
    enum real a75 = (-2187.0 / 6784.0);
    enum real a76 = (11.0 / 84.0);

    enum real c2 = (1.0 / 5.0);
    enum real c3 = (3.0 / 10.0);
    enum real c4 = (4.0 / 5.0);
    enum real c5 = (8.0 / 9.0);
    enum real c6 = (1.0);
    enum real c7 = (1.0);

    enum real b1 = (35.0 / 384.0);
    enum real b2 = (0.0);
    enum real b3 = (500.0 / 1113.0);
    enum real b4 = (125.0 / 192.0);
    enum real b5 = (-2187.0 / 6784.0);
    enum real b6 = (11.0 / 84.0);
    enum real b7 = (0.0);

    enum real b1p = (5179.0 / 57_600.0);
    enum real b2p = (0.0);
    enum real b3p = (7571.0 / 16_695.0);
    enum real b4p = (393.0 / 640.0);
    enum real b5p = (-92_097.0 / 339_200.0);
    enum real b6p = (187.0 / 2100.0);
    enum real b7p = (1.0 / 40.0);

    real[] K1, K2, K3, K4, K5, K6, K7;
    real x = x0;
    auto y = y0.dup;
    auto h = hmax;
    auto iter = maxiter;
    real[][] yout;
    yout ~= y0;
    real[] tout;
    tout ~= x0;

    while (x < x1)
    {
        /* Compute the function values */
        K1 = fxy(x, y, args);
        // writeln(K1, ' ', x);
        auto y1 = y.dup;
        y1[] += h * (a21 * K1[]);
        K2 = fxy(x + c2 * h, y1, args);
        auto y2 = y.dup;
        y2[] += h * (a31 * K1[] + a32 * K2[]);
        K3 = fxy(x + c3 * h, y2, args);
        auto y3 = y.dup;
        y3[] += h * (a41 * K1[] + a42 * K2[] + a43 * K3[]);
        K4 = fxy(x + c4 * h, y3, args);
        auto y4 = y.dup;
        y4[] += h * (a51 * K1[] + a52 * K2[] + a53 * K3[] + a54 * K4[]);
        K5 = fxy(x + c5 * h, y4, args);
        auto y5 = y.dup;
        y5[] += h * (a61 * K1[] + a62 * K2[] + a63 * K3[] + a64 * K4[] + a65 * K5[]);
        K6 = fxy(x + h, y5, args);
        auto y6 = y.dup;
        y6[] += h * (a71 * K1[] + a72 * K2[] + a73 * K3[] + a74 * K4[] + a75 * K5[] + a76 * K6[]);
        K7 = fxy(x + h, y6, args);

        auto error = y0.dup;
        error[] *= 0;
        error[] += (b1 - b1p) * K1[];
        error[] += (b3 - b3p) * K3[] + (b4 - b4p) * K4[];
        error[] += (b5 - b5p) * K5[] + (b6 - b6p) * K6[] + (b7 - b7p) * K7[];
        error = error.map!abs.array;

        // error control
        // auto delta = error.dup;
        // delta[] *= 0;
        // writeln(delta);
        // foreach (j, ref el; delta)
        // {
        //     el = 0.84 * pow(tol / (el+1e-16), 1 / 5.0);
        // }

        auto delta = 0.84 * pow(tol / (error.maxElement + 1e-16), (1.0 / 4.0));
        // writeln(error,h,delta);
        if (error.mean < tol)
        {
            x += h;
            y[] += h * (b1 * K1[] + b3 * K3[] + b4 * K4[] + b5 * K5[] + b6 * K6[]);
            writeln(y, h, K2, "\n");
            yout ~= y;
            tout ~= x;
        }

        if (delta <= 0.1)
        {
            h *= 0.1;
        }
        else if (delta >= 4.0)
        {
            h *= 4.0;
        }
        else
        {
            h *= delta;
        }

        if (h > hmax)
        {
            h = hmax;
        }

        if (x >= x1)
        {
            uint flag = 0;
            break;
        }

        else if (x + h > x1)
        {
            h = x1 - x;
        }

        else if (h < hmin)
        {
            uint flag = 1;
            break;
        }
        iter -= 1;
    }

    if (maxiter - iter <= 0)
    {
        uint flag = 2;
    }

    return tuple(tout, yout);
}

unittest
{
    real x0 = 0;
    real[] y0 = [1.0, 0.5, 2.0];
    real x1 = 100;
    real tol = 1.0e-5;
    real hmax = 1.0;
    real hmin = 0.01;
    uint maxiter = 1000;
    real[] fxy(real t, real[] Y, real[] P)
    {
        auto x = Y[0];
        auto y = Y[1];
        auto z = Y[2];
        real[] dy = [0, 0, 0];
        dy[0] += -x - y;
        dy[1] += x + P[0] * y;
        dy[2] += P[1] + z * (x - P[2]);
        return dy;
    }

    auto res = dopri(&fxy, x0, y0, x1, tol, hmax, hmin, maxiter, [0.1, 0.1, 14]);
    // writeln(res[1]);
}
