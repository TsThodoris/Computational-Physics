#---------- Monty Hall Problem ----------#
# Monty Hall problem for three cases:
#Case 1: stick with initial choice
#Case 2: always change door
#Case 3: always choose door A, if C opens switch to B , else nothing
import numpy,random 
from numpy.random import default_rng
import matplotlib.pyplot as plt

#----- Doors
doors = ["A","B","C"]

#----- Parameters and arrays

ns = numpy.linspace(10,1000000,100)   #iterations
gen = numpy.random.default_rng(12345)   #random number generator
case_val = 0    #values: 0,0.5 and 1 for case a,b and c respectively
data = numpy.zeros(len(ns))
plt.figure(1)

def set_door(x):
    if x<=1/3: y = "A"
    elif (x>1/3 and x<=2/3): y = "B"
    elif (x>2/3 and x<=3/3): y = "C"

    return y

for ind,n in enumerate(ns):
    c1 = 0;
    c2 = 0;
    c3 = 0;
    for i in range(round(n)):
    
        if case_val == 0: second_chance = 0
        elif case_val == 0.5: second_chance = 1
    
        prize_num = gen.random()
        player_num = gen.random()
        #strings
        prize = set_door(prize_num)
        selection = set_door(player_num)
        if (case_val == 1 and selection!="A"): selection = "A"

        if selection == prize:
            res = list(set(doors)-set(prize))
            monty = random.choice(list(set(doors)-set(random.choice(res)))) #monty opens one of the doors without prize
            alternative = random.choice(list(set(doors)-set(monty)-set(prize)))
        else:
            monty = random.choice(list(set(doors)-set(selection)-set(prize)))
            alternative = random.choice(list(set(doors)-set(monty)-set(selection)))

   
        if (case_val == 1 and monty == "C"): second_chance = 0.5
        elif(case_val==1 and monty!="C"): second_chance = 0
 
        if second_chance == 1:
            if alternative == prize: c2+=1
        elif second_chance == 0:
            if (selection == prize and case_val == 0): c1+=1
            elif (selection == prize and case_val == 1): c3+=1
        elif second_chance ==0.5:
            if alternative == prize: c3+=1
      
    if case_val ==0: c = c1
    elif case_val == 0.5: c = c2
    elif case_val == 1: c = c3
    data[ind] = c/n

plt.plot(ns,data)
plt.scatter(ns,data)
plt.title("Probability of wining with increasing number of iterations")
plt.xlabel("n (iterations)")
plt.ylabel("Probability")
plt.figure(2)
plt.hist(data)
plt.xlabel("Probability")
plt.xlim([0.3,0.4])
plt.title("Histogram of probability values for each number of iterations")
plt.show()
