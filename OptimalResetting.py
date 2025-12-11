'''
This codes simulates free Brownian motion under stochastic resetting via optimal retrun protocol.
It is structured as loops studying several choices of parameters : various stiffness of the return potential
and various protocol duration
The physical parameters chosen correspond to typical colloidal experiments [https://doi.org/10.48550/arXiv.2306.09503]
It outputs the measured values of work and mean first passage time, as well as their statistical error
The associated plot.m code ensures the visual rendering of these numerical resutls
'''
import numpy as np
import matplotlib
import matplotlib.pyplot as plt
import pandas as pd
import trackpy as tp
from tqdm import tqdm
from scipy.io import savemat
import os

kB = 1.3806e-23 # Boltzmann constant
T = 300 # room temperature
gamma = 2.8e-8 # Stokes viscous drag
D = kB*T/gamma # Diffusion coefficient
l = 0.1e-6 # Location of the targer
r = 36  # resetting rate
dt = 1e-5 # time-step
Ntraj = 300 # number of trajectories
NCaseTf = 7 # number of cases of protocol length
KappaValue = [1, 5, 25, 50, 100] # in pN/\mu m
NCaseKappa = np.size(KappaValue)
MinTimeArray = np.linspace(0.005, 0.5, NCaseKappa)
MaxTimeArray = np.linspace(5, 35, NCaseKappa)

def lambdaOpt(t, tf, ui, k, gamma, lambdaf):
    # optimal translation protocol starting at the position ui (here, measured position)
    # and ending up at lambdaf (here, prescribed destination : 0), within a finite time tf
    return(ui + (k*t/gamma + 1) * (lambdaf-ui)/(k*tf/gamma))

for ll in range(NCaseKappa): # loop over the various stiffness [KappaValue]
    print(KappaValue[ll])
    k = KappaValue[ll] * 1e-6
    Alltf = gamma / k * np.logspace(np.log10(MinTimeArray[ll]), np.log10(MaxTimeArray[ll]), NCaseTf)
    AllMeanWork = []
    AllErrWork = []
    AllMeanfpt = []
    AllMeanfptInst = []
    AllErrfpt = []
    for ii in range(NCaseTf): # loop over the protocol durations [Alltf]
        print(str(ii+1) + ' / ' + str(NCaseTf))
        tf = Alltf[ii]
        AllWork = []
        AllFpt = []
        AllFptInst = []
        for i in tqdm(range(Ntraj)): # loop over trajectories
            traj = []
            protocol = []
            x = 0
            t = 0
            count = 0
            work = 0
            works = [] 
            while x < l: # until the particle reaches the targer
                q = np.random.binomial(1, 1-np.exp(-r * dt)) # binomial random variable triggering reset
                if q == 1: # if resetting happens: reset protocol
                    count += 1
                    x0 = np.copy(x) # position when resetting starts
                    work +=  0.5 * k * ( x0 -  lambdaOpt(0, tf, x0, k, gamma, 0) )**2 # work cost of the first discontinuity
                    for j in range(int(tf/dt)):
                        t += dt
                        lam = lambdaOpt(t, tf, x0, k, gamma, 0) # optimal return protocol imposed
                        protocol.append(lam)
                        x += - k/gamma * (x - lam) * dt + np.sqrt(2 * kB * T * dt / gamma) * np.random.normal(0,1) # evolution during return
                        traj.append(x)
                        work += -  x0 / tf * k *(lam - x) * dt # work cost of the return protocol
                    work += 0.5 * k * x**2 - 0.5 * k *( x -  lambdaOpt(tf, tf, x0, k, gamma, 0) )**2 # work cost of the final discontinuity
                    if work != 0:
                        works.append(work)
                    t = 0
                    work = 0 
                else: # if resetting does not happen : free diffusion
                    protocol.append(0)
                    x += np.sqrt(2 * kB * T * dt / gamma) * np.random.normal(0,1) # free diffusion
                    traj.append(x)

            if len(works) != 0:
                for wval in works:
                    AllWork.append(wval)
            AllFpt.append(len(traj)*dt) # first passage time including the finite time returns
            AllFptInst.append(len(traj)*dt - count*tf) # measured ideal instantaneous first passage time
        MeanWork = 0
        wcnt = 0
        for ww in AllWork:
            if ww != 'nan':
                MeanWork += ww
                wcnt += 1
        MeanWork = MeanWork/wcnt # mean work cost of the trajectory, until first passage

        ErrWork = np.sqrt(np.var(AllWork))/np.sqrt(Ntraj)
        MeanFpt = np.mean(AllFpt)
        ErrFpt = np.sqrt(np.var(AllFpt))/np.sqrt(Ntraj)
        MeanFptInst = np.mean(AllFptInst)
        ErrFptInst = np.sqrt(np.var(AllFptInst))/np.sqrt(Ntraj)
        AllMeanWork.append(MeanWork)
        AllMeanfpt.append(MeanFpt)
        AllMeanfptInst.append(MeanFptInst)
        AllErrWork.append(ErrWork)
        AllErrfpt.append(ErrFpt)

    z = l * np.sqrt(r / D)
    fudge = 1 - 1 / (1 - np.exp(z)) + z**2 / (4 * (1 - np.cosh(z * np.sqrt(1 - np.exp(-z))))) # correction to the second moment due to absorbing boundary condition in l
    savemat("AllMeanWork" + str(KappaValue[ll]) + ".mat", {'AllMeanWork': AllMeanWork})
    savemat('Fuge.mat', {'fudge': fudge})
    savemat("Alltf" + str(KappaValue[ll]) + ".mat", {'Alltf': Alltf})
    savemat("AllMeanfpt" + str(KappaValue[ll]) + ".mat", {'AllMeanfpt': AllMeanfpt})
    savemat("AllMeanfptInst" + str(KappaValue[ll]) + ".mat", {'AllMeanfptInst': AllMeanfptInst})
    savemat("AllErrWork" + str(KappaValue[ll]) + ".mat", {'AllErrWork': AllErrWork})
    savemat("AllErrfpt" + str(KappaValue[ll]) + ".mat", {'AllErrfpt': AllErrfpt})
