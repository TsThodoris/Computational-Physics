#-------------------- Project Quantum Mechanics
#-------------------- Tsitlakidis Theodoros 4409

import numpy as np, scipy as sc, matplotlib.pyplot as plt
from math import *
from scipy.integrate import solve_ivp
from scipy.integrate import simpson

#First of all we write the Lane-Emden equation as a system of two first order
def le(y,x,g):
    u,z = x
    return [z,-(2/y)*z-(u**(1/(g-1)))]

#---------- Parameters
tmin = 0.00001
tmax = 30
t_span = np.array([tmin,tmax])#t_span
step = 0.001     #time step
gamma_mat = np.linspace(1.2,1.9,169)
#---- Parameters
hbar = 197.327  #hbarc [MeV fm]
m_N = 938.27  #nucleon mass
m_e = 0.511   #electron mass
mu = 2   #num on nucleons per electron
G = 3.743*10**(19) #in MeV
n0 = 0.16144 #nuclear density [fm^-3]
K1 = ((hbar**2/(15*m_e*pi**2))*(3*pi**2/(m_N*mu)))**(5/3)   #small mass
K2 = ((hbar/(12*pi**2))*(3*pi**2/(m_N*mu)))**(4/3)      #large mass
M_rescale = 200
S = np.zeros(len(gamma_mat))
S3 = np.zeros((3,len(gamma_mat)))
M = np.zeros(len(gamma_mat))
ic = np.array([1,0])
plot_le = np.linspace(1.2,1.9,8)
kappa_val = np.array([0.95,1.00,1.05])
iterations = 10000


# Figure 2 Normalized modal fraction f(|k|) for sample values of polytropic index gamma (figure 1 of paper)
plt.figure(2)
plt.ylabel(r"$\bar f(|k|)$")
plt.xlabel(r"$k/\sqrt{4\pi G/K\rho_0^{\gamma-2}}$")
plt.xlim([0, 1.5])
plt.ylim([0, 1.1])  
plt.grid()

# Figure 4 Configurational entropy times \rho^(-1) (continuous line) and mass (dotted line) versus polytropic 
# index gamma (figure 2 of paper)
plt.figure(4)
plt.xlabel(r"$\gamma$")
plt.xlim([1.25, 1.7])
plt.ylim([0.4, 1.3])
plt.grid()

# Figure 3 Configurational entropy versus polytropic index gamma for polytropes. We display results for several 
# choices of cutoff for k_min (figure 5 of paper)
plt.figure(3)
plt.xlabel(r"$\gamma$")
plt.ylabel(r"$Sa^{3}$")
plt.xlim([1.25, 1.75])
plt.ylim([4.6, 5.6])
plt.axvline(x=4/3,color='k',ls='--')
plt.text(4/3, 4.63, "4/3", rotation=0)
plt.axvline(x=5/3,color='k',ls='--')
plt.text(5/3, 4.63, "5/3", rotation=0)
plt.grid()

#---------- Solving 
for i,gamma in enumerate(gamma_mat):   
    eval_i = np.linspace(tmin,tmax,iterations)
    sol = solve_ivp(le,t_span,ic,method='RK45',args=(gamma,),t_eval = eval_i)
    X = sol.t
    Y =sol.y
    pos_theta = np.where(Y[0]>=0)
    thetas = Y[0][pos_theta]
    ksi = X[pos_theta]
    
    if gamma in plot_le:
        num = str(round(gamma,1))
        l = 'γ='+num
        plt.figure(1)
        plt.plot(X,Y[0],label=l)
        plt.legend()
#---------- Plotting
        plt.ylim([0,1])
        plt.xlabel('ξ')
        plt.ylabel('θ\'(ξ)')
        plt.title('Lane-Emden solution for various γ')
#---------- Continuing with every else
    min_k = (pi/ksi[-1])/kappa_val
    #computing of h(kmin) and h(k) 
    for inde,j1 in enumerate(min_k):
        k = np.linspace(j1,100*j1,iterations)
        func = thetas**(1/(gamma-1))*np.sin(j1*ksi)*ksi
        hkmin = (simpson(func,ksi,dx=0.001)*(1/j1))**2
        hk_mat = np.zeros(iterations)
        for index,jj in enumerate(k):
            func_jj = thetas**(1/(gamma-1))*np.sin(jj*ksi)*ksi
            hk_mat[index] = (simpson(func_jj,ksi,dx=0.001)*(1/jj))**2
        f = hk_mat/hkmin

        mass = (thetas**(1/(gamma-1)))*ksi**2
        M[i] = (4*pi*((gamma/(gamma-1))**((3/2)))*simpson(mass,ksi,dx=0.001))/M_rescale

        #Configuration entropy computing
        s_func = f*np.log(f)*k**2
        if inde == 1:
            S[i] = -4*pi*((gamma/(gamma-1)))**(-3/2)*simpson(s_func,k,dx=0.0001)
        S3[inde][i] = -4*pi*simpson(s_func,k,dx=0.0001)


        if inde == 1:
            plt.figure(2)
            index, = np.where(f<=1)
            if gamma==1.2:
                plt.scatter(k[index[0]]/np.sqrt(gamma/(gamma-1)),f[index[0]],marker='v',color='r')
                plt.plot(k[index[0]:]/np.sqrt(gamma/(gamma-1)),f[index[0]:],'--',color='r',label=f'γ = {gamma}')
            elif gamma==1.4:
                plt.scatter(k[index[0]]/np.sqrt(gamma/(gamma-1)),f[index[0]],marker='v',color='g')
                plt.plot(k[index[0]:]/np.sqrt(gamma/(gamma-1)),f[index[0]:],color='g',label=f'γ = {gamma}')
            elif gamma==1.7:
                plt.scatter(k[index[0]]/np.sqrt(gamma/(gamma-1)),f[index[0]],marker='v',color='c')
                plt.plot(k[index[0]:]/np.sqrt(gamma/(gamma-1)),f[index[0]:],'-.',color='b',label=f'γ = {gamma}')
    
                                
plt.figure(2)
plt.title('f(|k|) - k')
plt.legend

plt.figure(3)
plt.title('$Sa^{3}$ - γ')
plt.scatter(gamma_mat[np.argmax(S3[0])],np.amax(S3[0]),marker='v',color='r')
plt.scatter(gamma_mat[np.argmax(S3[1])],np.amax(S3[1]),marker='v',color='r')
plt.scatter(gamma_mat[np.argmax(S3[2])],np.amax(S3[2]),marker='v',color='r')
plt.scatter(gamma_mat[np.argmin(S3[0])],np.amin(S3[0]),marker='o',color='b')
plt.scatter(gamma_mat[np.argmin(S3[1])],np.amin(S3[1]),marker='o',color='b')
plt.scatter(gamma_mat[np.argmin(S3[2])],np.amin(S3[2]),marker='o',color='b')
plt.plot(gamma_mat,S3[0],'--',color='b',label=r'$\pi/(0.95R)$')
plt.plot(gamma_mat,S3[1],color='g',label=r'$\pi/(1.00R)$')
plt.plot(gamma_mat,S3[2],'-.',color='r',label=r'$\pi/(1.05R)$')
plt.legend()


plt.figure(4)
plt.title('S and M - γ')
plt.plot(gamma_mat,S,label=r'$S\rho_0^{-1}/\left( \left( \frac{K}{4\pi G}\right)^{-\frac{3}{2}} \rho_c^{2-\frac{3}{2}\gamma} \right)$')
plt.plot(gamma_mat,M,'--',color='g',label=r'$M/\left( 200\left( \frac{K}{4\pi G}\right)^{\frac{3}{2}} \rho_c^{\frac{3}{2}\gamma -2} \right)$')
plt.legend()


plt.show()


