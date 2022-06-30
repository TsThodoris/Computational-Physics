#---------- 1D Ising Model with Metropolis Algorithm ----------#

import numpy as np
import random
from numpy.random import default_rng
import matplotlib.pyplot as plt


gen = np.random.default_rng(12345)  #random number generator

def initial_state(nsteps):
    ini_state = 2*gen.integers(0,2,nsteps)-1
    return ini_state

def calc_energy(state,J):
    counter = len(state)
    E = 0.0
    for i in range(counter):
        si = state[i]
        sj = state[(i+1)%counter]   #% for indexing
        E = E-J*(si*sj)
    return E

def energy_diff(J,S,sl,sr):
    return 2*J*S*(sl+sr)
    
def metropolis(nsteps,nsites,beta,J):
    state = initial_state(nsteps)
    E = calc_energy(state,J)
    spins = np.zeros(nsteps)
    energy_mat = np.zeros(nsteps)
    
    for i in range(nsteps):
        ind = gen.integers(0,nsteps)
        change_spin = state[ind]
        si = state[(ind-1)%nsites]
        sj = state[(ind+1)%nsites]

        dE = energy_diff(J,change_spin,si,sj)
        r = gen.uniform()
        #Metropolis
        if r < np.exp(-dE*beta):
            state[ind] *= -1    #accept change
            E += dE
        spins[i] = state.mean()
        energy_mat[i] = E

    return spins,energy_mat

def plotting(func,fig_num,tl,yl,xl):
    axs[fig_num].plot(func)
    axs[fig_num].set_title(tl)
    
    

#Sampling
nsteps = 100000
nsites = 10000  #lattice sites
J = 1  # stregnth of exchange interaction
T_mat = np.array([1,10,20])
fig, axs =plt.subplots(3,sharex = True)
fig.tight_layout()
for ind,T in enumerate(T_mat):
    beta = 1/T
    spins,energies = metropolis(nsteps,nsites,beta,J)
    title = r'$\beta=%.3f,T=%d$'%(beta,T)
    xls = 'Steps'
    yls = r'$m$'
    plotting(spins,ind,title,yls,xls)
plt.xlabel('Steps')
fig.supylabel(r'$m$')
plt.show()
