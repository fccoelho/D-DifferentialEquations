# Differential Equations Solver Suite for D 

This repo aims to be an open initiative to develop a suite of Solvers fully implemented the [D language](https://dlang.org). This is by no means original work, and should consist in porting and adapting existing implementation in other languages. Interesting sources of reference implementations are:
 - [Scipy's IVP Python implementation](https://github.com/scipy/scipy/tree/master/scipy/integrate/_ivp). Should be easy to port.
 - Add your source...

The current code is just a proof of concept, to be completely replaced by better implementations.

## Getting started
This library is in its early stages of development, so you can clone the repo, and run the tests:
```bash
$ dub test -f
Generating test runner configuration 'drunge-kutta-test-library' for 'library' (library).
Performing "unittest" build using /usr/bin/dmd for x86_64.
drunge-kutta ~master: building configuration "drunge-kutta-test-library"...
Linking...
Running ./drunge-kutta-test-library 
 ✓ solvers rk4 basic test
 ✓ solvers Dopri basic test
```
Or, to sun some early benchmarks: 
```bash
dub
```

You can also validate the implemented models against Scipy's `odeint`. This script calculates the mean-squared-error between the two solutions:

```bash
$ python validate.py
SIR MSE:  1.3638751214735066e-05
```
![SIR](https://github.com/fccoelho/D-DifferentialEquations/blob/master/validate_sir.png)


## Contributing
I have no intentions of moving this library forward on my own, so any help is very much welcome! just open an issue with some topic you would like to contribute, fork the repo and start to work!

We have also enabled [discussions](https://github.com/fccoelho/D-DifferentialEquations/discussions) here, have a try!
