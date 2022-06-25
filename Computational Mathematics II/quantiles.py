#---------- Quantile Plot ----------#

import numpy,random
from numpy.random import default_rng
import matplotlib.pyplot as plt
import statsmodels.api as sm
import scipy.stats as stats


n = 500
gen = numpy.random.default_rng(12345)
rand_nums = numpy.zeros(n)
inv_nums = numpy.zeros(n)
lamda = 0.5

for i in range(n):
    rand_nums[i] = gen.exponential(scale=1.0)

    inv_nums[i] = -numpy.log(1-rand_nums[i])/lamda
    
inv_nums = inv_nums[numpy.logical_not(numpy.isnan(inv_nums))]
#ploting
sm.qqplot(rand_nums,line='45',fit=True,dist=stats.norm)
plt.title("Normal Distribution and exp random numbers")
sm.qqplot(inv_nums,line='45',fit=True,dist=stats.expon)
plt.title("Exponential Distribution and inverse exp random numbers")
sm.qqplot(inv_nums,line='45',fit=True,dist=stats.norm)
plt.title("Normal Distribution and inverse exp random numbers")
sm.qqplot(inv_nums,line='45',fit=True,dist=stats.uniform)
plt.title("Uniform Distribution and invers exp random numbers")
plt.show()
