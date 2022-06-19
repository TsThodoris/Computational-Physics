import numpy as np
import math as ma
import matplotlib.pyplot as plt
import matplotlib.patches as mpatches

p0=20
m0=0
r0=10**(-5)
e0=64.5787*p0**(0.437371)+64.5754*p0**(0.427547)               #######Changes for different EOS#######


def dpdr(E,M,P,r):
    return -1.474*((E*M)/(r**2))*(1+(P/E))*(1+11.2*(10**(-6))*(r**3)*(P/M))*((1-2.948*(M/r))**(-1))

def dmdr(E,r):
    return 11.2*10**(-6)*r**2*E

######EOS
def eos1(P):       ###mx=1.0mn
    return 14.5809*P**(0.813358)+88.4385*P**(0.42247)

def eos2(P):        ###mx=1.2mn
    if P>=0.0001 and P<=45.9396:
        return 0.625425*P**(1.02671)+92.825*P**(0.34556)
    else:
       return 18.7993*P**(0.746537)+18.7994*P**(0.746537)

def eos3(P):        ###mx=1.5mn
    if P>=0.0001 and P<=190.008:
        return 19.7577*P**(0.565022)+24.5308*P**(0.275924)+24.5353*P**(0.272785)+24.5356*P**(0.270283)
    else:
        return 7.74685*P**(0.889116)+7.74685*P**(0.889116)
    

def eos4(P):        ###mx=2.0mn
    return 20.6993*P**(0.791071)+20.6993*P**(0.791071)
   

def eos_5(P):       ###MDI-1
    if P>=0.0001 and P<=559.174:
        return 0.879941*P**(0.953308)+92.4507*P**(0.344931)
    else:
        return 4.1844*P**(0.81449)+95.00135*P**(0.31736)

def eos_6(P):       ###HHJ-1
    return 1.78429*P**(0.93761)+106.93652*P**(0.31715)

def eos_7(P):       ###NLD
    return 119.05736+304.80445*(1-np.exp(-P/48.61465))+33722.34448*(1-np.exp(-P/17499.47411))

def eos_dminter1(P):
    return 195.823*P**(0.287184)+2.65987*P**(0.820264)

def eos_dminter2(P):
    return 0.120481*P**(1.21764)+273.762*P**(0.235691)


def eos_dminter3(P):    #z=50
    return 2.10765*P**(0.837683)+162.348*P**(0.311779)

def eos_dminter4(P):    #z=100
    return 67.2235*P**(0.409544)+67.2249*P**(0.383884)

def eos_dminter5(P):
    return 64.8854*P**(0.455944)+64.58854*P**(0.449277)

def eos_dmintter6(P):
    return 62.1396*P**(0.424558)+62.1397*P**(0.416542)

def eos_dminter7(P):
   return 51.67*P**(0.496417)+51.6706*P**(0.489301)

def eos_dminter8(P):
    return 64.5787*P**(0.437371)+64.5754*P**(0.427547)

#### Function for solving
def eos(P):
    c0=31.93753
    c1=10.82611*np.log10(P)
    c2=1.29312*(np.log10(P)**2)
    c3=0.08014*(np.log10(P)**3)
    c4=0.00242*(np.log10(P)**4)
    c5=0.000028*(np.log10(P)**5)
    if P>9.34975*10**(-5) and P<=0.184:
        return 0.00873+103.17338*(1-ma.exp(-P/0.38527))+7.34979*(1-ma.exp(-P/0.01211))
    elif P>4.1725*(10**(-8)) and P<=9.34975*(10**(-5)):
        return 0.00015+0.00203*(1-ma.exp(-P*344827.5))+0.10851*(1-ma.exp(-P*7692.3076))
    elif P>1.44875*(10**(-11)) and P<=4.1725*(10**(-8)):
        return 0.0000051*(1-ma.exp(-P*0.2373*(10**(10))))+0.00014*(1-ma.exp(-P*0.4020*(10**8)))
    elif P<=1.44875*(10**(-11)):
        return 10**(c0+c1+c2+c3+c4+c5)
    else:
        return 64.5787*P**(0.437371)+64.5754*P**(0.427547)                        ######Changes for different EOS########


h=0.001          ###h=0.001
Mf=np.array([])                 ####Final Mass,Pressure and Radius
Pf=np.array([])
Rf=np.array([])
radius=np.array([r0])
mass=np.array([m0])
press=np.array([p0])
enr=np.array([e0])
m1=11.2*10**(-6)*h**3*e0
mass=np.append(mass,m1)
    
for i in range(5,1200,5):
  p0=i
  m0=0
  r0=10**(-5)
  e0=64.5787*p0**(0.437371)+64.5754*p0**(0.427547)                                ######Changes for different EOS######
  radius=np.array([r0])
  mass=np.array([m0])
  press=np.array([p0])
  enr=np.array([e0])
  m1=11.2*10**(-6)*h**3*e0
  mass=np.append(mass,m1)
  for j in range(1,30000):                                                  #####The number of iterations may have to increase for different EOS#####
    k1=dpdr(e0,m1,p0,r0)
    l1=dmdr(e0,r0)
    k2=dpdr(e0+k1*h/2,m1+l1*h/2,p0+k1*h/2,r0+h/2)
    l2=dmdr(e0+k2*h/2,r0+h/2)
    k3=dpdr(e0+k2*h/2,m1+l2*h/2,p0+k2*h/2,r0+h/2)
    l3=dmdr(e0+k3*h/2,r0+h/2)
    k4=dpdr(e0+k3*h,m1+l3*h,p0+k3*h,r0+h)
    l4=dmdr(e0+k3*h,r0+h)
    p_new=p0+(k1+2*k2+2*k3+k4)*h/6
    m_new=m1+(l1+2*l2+2*l3+l4)*h/6
    e0=eos(p0)
    r0=r0+h
    p0=p_new
    m1=m_new
    mass=np.append(mass,m_new)
    press=np.append(press,p_new)
    radius=np.append(radius,r0)
    enr=np.append(enr,e0)
  j=0
  while press[j]>0:
    j=j+1
  radii=radius[j].real
  Rf=np.append(Rf,radii)
  mf=mass[j].real
  Mf=np.append(Mf,mf)
  
  
#####Ploting
  
plot_patch=mpatches.Patch(color='orange',label='$7.74685P^{0.889116}+7.74685P^{0.889116}$')

plt.xlabel("R [Km]")
plt.ylabel("M [$M_{\odot}]$")
plt.legend(handles=[plot_patch])
plt.plot(Rf,Mf,'orange')
plt.show()
print(Rf,Mf)


    














