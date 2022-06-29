#---------- Accept/Reject Algorithm ----------#

import numpy,random,math
from numpy.random import default_rng
import matplotlib.pyplot as plt

n = 100
gen = numpy.random.default_rng(12345)

#The examples disttribution
def f(x):
    return (2/numpy.pi)*((1+x)**(-1/2))*((1-x)**(-1/2))

#plot of function f(x)
plt.figure(1)
xs = numpy.linspace(0,1,100)
ys = f(xs)
plt.plot(xs,ys,'orange',label=r'f(x)=$\frac{2}{\pi}\frac{1}{(1+x)^{1/2}(1-x)^{1/2}}$')
ys = [i for i in ys if not math.isinf(i)]
f_max = max(ys)

  
def accept_reject(n):
    accept = numpy.zeros(n)
    u = 0
    while u < n:
        t = gen.uniform(0,1)
        y = gen.uniform(0,f_max)
        if y < f(t):
            accept[u] = t
            u+=1
    #remove zero values
    accept = [k for k in accept if k!=0]
    return accept

a = accept_reject(100*n)
plt.hist(a,bins=100,density=True,label='Sample Distribution')
plt.title("Acceptance Plot")
plt.xlabel('x')
plt.ylabel('y')
plt.legend()


#Plot of total sampling and mask of accepted values
t = gen.uniform(0,1,100*n)
y = gen.uniform(0,f_max,100*n)
plt.figure(2)
plt.scatter(t,y,s=1,c='blue',label='total sampling')
mask = y<f(t)
plt.scatter(t[mask],y[mask],s=1.2,c='orange',label='mask')
plt.title('Acceptance Plot with mask of total sampling')
plt.xlabel('x')
plt.ylabel('y')
plt.show()

    
        
