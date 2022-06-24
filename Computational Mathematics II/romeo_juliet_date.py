#---------- Romeo and Juliet Date ----------#
import numpy,random
from numpy.random import default_rng
from statistics import *
import matplotlib.pyplot as plt

gen = numpy.random.default_rng(12345)
plt.figure(1)
ns = numpy.arange(10,1000000,1000)
meet_data = numpy.zeros(len(ns))
romeo = numpy.zeros(len(ns))
juliet = numpy.zeros(len(ns))
same_time_mat = numpy.zeros(len(ns))

for ind,n in enumerate(ns):
    r_num = 0
    j_num = 0
    same_time = 0;
    meet = 0;
    for i in range(n):
        romeo_arrival = gen.random()
        juliet_arrival = gen.random()
        #random number correspondace to time of arrival
        rt = 60*romeo_arrival
        jt = 60*juliet_arrival
        
        if rt<jt:
            r_num+=1
            if rt+15>jt: meet+=1
        elif rt>jt:
            j_num+=1
            if jt+15>rt: meet+=1
        else:
            same_time+=1
            meet+=1
    meet_percentage = meet/n
    romeo_first = r_num/n
    juliet_first = j_num/n
    meet_data[ind] = meet_percentage
    romeo[ind] = romeo_first
    juliet[ind] = juliet_first
    same_time_mat[ind] = same_time/n

plt.plot(ns,meet_data,label="data")
plt.plot([0,max(ns)],[7/16,7/16],label="7/16")
plt.xlabel("n (iterations)")
plt.ylabel("Probability of meeting")
plt.title("Romeo and Juliet Date Probability")
plt.legend()
plt.figure(2)
plt.hist(meet_data)
plt.plot([7/16,7/16],[0,100],label="7/16")
plt.title("Probability of meeting")
plt.xlabel("Probability of Meeting")
plt.legend()
plt.show()

print("Probability of Romeo arriving first:",mean(romeo))
print("\n","Probability of Juliet arriving first:",mean(juliet))
print("\n","Probability of arriving at the same time:",mean(same_time_mat))
    
        


