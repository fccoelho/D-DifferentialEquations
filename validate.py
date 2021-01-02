"""
This script checks the results of the benchmark models in `App.d` by comparing them
to the results obtained with Scipy Odeint
"""

from scipy.integrate import odeint
import os
import numpy as np
import pylab as P

## Generates the D solutions
os.system('dub')

## Python models

def sir(Y, t, *pars):
    s,i,r = Y
    beta, gamma = pars
    return [
        -beta*i*s, 
        beta*i*s-gamma*i, 
        gamma*i
    ]
trange = np.arange(0,500,0.01)
res_sir = odeint(sir,[0.999,0.001,0],trange,(0.3,0.05))

sir_d = np.genfromtxt('sir.csv',delimiter=',',skip_header=True);
print("SIR MSE: ",((sir_d[:,1:]-res_sir)**2).mean())
P.plot(sir_d[:,0], sir_d[:,2], '-*',label="D")
P.plot(trange, res_sir[:,1],label="Python")
P.title("SIR model")
P.grid()
P.legend()
P.savefig("validate_sir.png")
P.show()
