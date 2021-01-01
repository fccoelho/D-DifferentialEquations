import std.stdio;
import std.array;
import std.algorithm;
import std.conv;
import std.file;
import std.datetime.stopwatch;
import core.exception : RangeError;

import solvers; // Our solvers

/**
Save the simulation as a CSV file
*/
void save(T)(string filename, string[] header, T data)
{
    File outf = File(filename, "w");
    auto head = join(header, ",");
    outf.writeln(head);
    foreach (i, real[] row; data[1])
    {
        try
        {
            outf.writeln(text(data[0][i]), ",", row.map!text.join(","));
        }
        catch (RangeError e)
        {
            writefln("some error: ", row);
        }
    }
    outf.close();
}

/**
Examples
*/

/**
SIR model
*/
void sir(string method="rk4"){
       // model function
    writeln("SIR model.");
    real[] model(real t, real[] Y, real[] P)
    {
        auto s = Y[0];
        auto i = Y[1];
        auto r = Y[2];
        const real beta = P[0];
        const real gamma = P[1];

        return [-beta*s*i, beta*s*i - (gamma*i), gamma*i];
    }

    real[] inits = [1.0, 0.001, 0.0];
    real[] times;
    times ~= 0;
    foreach (i; 0 .. 10_000)
    {
        times ~= times[$ - 1] + 0.01;
    }
    auto sw = StopWatch(AutoStart.no);
    sw.start();
    auto res = odeint(&model, inits, times, [0.3, 0.05], method,1e-6,0.01);
    sw.stop();
    writefln("Time of the Model run: %s seconds", sw.peek());
    save("sir.csv", ["t", "S", "I", "R"], res);
}
unittest{

}
/**example 1*/
void example_1()
{
    writeln("Example dodeint run.");
    real[] fun(real t, real[] Y, real[] P)
    {
        const auto u = Y[0];
        const auto v = Y[1];
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
    auto sw = StopWatch(AutoStart.no);
    sw.start();
    auto res = rk4(&fun, inits, times, 1e-6, [0.1, 2.0, 0.2]);
    sw.stop();
    writefln("Time of the Mode run: %s seconds", sw.peek());
    save("res.csv", ["u", "v"], res);
}

///Rossler function
void rossler(string method = "rk4")
{
    // model function
    writeln("Rossler oscillator.");
    real[] osc(real t, real[] Y, real[] P)
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

    real[] inits = [1.0, 0.5, 2.0];
    real[] times;
    times ~= 0;
    foreach (i; 0 .. 3333)
    {
        times ~= times[$ - 1] + 0.15;
    }
    auto sw = StopWatch(AutoStart.no);
    sw.start();
    auto res = odeint(&osc, inits, times, [0.1, 0.1, 14], method);
    sw.stop();
    writefln("Time of the Mode run: %s seconds", sw.peek());
    save("rossler.csv", ["t", "x", "y", "z"], res);
}

void main()
{

    rossler("rk4");
    sir("dopri");
    // example_1();
    // auto dmd = execute(["graph", "-T", "res.csv"]);
}
