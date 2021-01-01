import std.stdio;
import solvers;
import std.array;
import std.algorithm;
import std.conv;
import std.file;
import std.datetime.stopwatch;
import core.exception : RangeError;

void save(T)(string filename, string[] header, T data)
{
    File outf = File(filename, "w");
    auto head = join(header, ",");
    outf.writeln(head);
    foreach (i, real[] row; data)
    {
        try
        {
            outf.writeln(row.map!text.join(","));
        }
        catch (RangeError e)
        {
            writefln("some error: ", row);
        }
    }
    outf.close();
}

void example_1()
{
    writeln("Example dodeint run.");
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
    auto sw = StopWatch(AutoStart.no);
    sw.start();
    auto res = dodeint(&fun, inits, times, 1e-6, [0.1, 2.0, 0.2]);
    sw.stop();
    writefln("Time of the Mode run: %s secconds", sw.peek());
    save("res.csv", ["u", "v"], res);
}

///Rossler function
void rossler(){
    // model function
    writeln("Rossler oscillator.");
    real[] fun(real t, real[] Y, real[] P)
    {
        auto x = Y[0];
        auto y = Y[1];
        auto z = Y[2];
        real[] dy = [0, 0, 0];
        dy[0] += -x-y;
        dy[1] += x + P[0] * y;
        dy[2] += P[1] + z*(x-P[2]);
        return dy;
    }
    real[] inits = [1.0,0.5, 2.0];
    real[] times;
    times ~= 0;
    foreach (i; 0 .. 3333)
    {
        times ~= times[$ - 1] + 0.15;
    }
    auto sw = StopWatch(AutoStart.no);
    sw.start();
    auto res = dodeint(&fun, inits, times, 1e-6, [0.1, 0.1, 14]);
    sw.stop();
    writefln("Time of the Mode run: %s secconds", sw.peek());
    save("rossler.csv", ["x", "y","z"], res);
}
void main()
{
    
    rossler();
    // example_1();
    // auto dmd = execute(["graph", "-T", "res.csv"]);
}
