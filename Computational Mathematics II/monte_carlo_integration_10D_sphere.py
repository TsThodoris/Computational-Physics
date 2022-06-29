#---------- Monte Carlo Integration for n-D Sphere Volume ----------#

import numpy,random
import numpy.linalg as lin
from numpy.random import default_rng

n = 1000000    #iterations
dim = 10   #dimensions
gen = numpy.random.default_rng(12345)   #generator of random numbers

#Cube of side (i.e.) 2
s = 2
nums = gen.uniform(-s/2,s/2,(n,dim))
vol = s**dim
R = lin.norm(nums,axis=1)
u = numpy.sum(R<=s/2)
Volume = (vol/n)*u

print("The Volume of the ",dim,"D sphere, of radius ",s,"is:",Volume)
