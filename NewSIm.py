
import numpy as np
import math 
from scipy.integrate import odeint 
import random
import matplotlib.pyplot as plt 


newdir = "C:/Users/Daniel Park/Desktop/SimGraphs/"

for i in range(100):
    umax = 2
    dc = random.uniform(0, random.uniform(0,2))
    ks = random.uniform(0, random.uniform(0,2))
    ki = random.uniform(0, random.uniform(0,2))
    cmax = 0.225
    g = random.uniform(0, random.uniform(0,2))
    dp = random.uniform(0, random.uniform(0,2))
    w = 0.06
    qmax = random.uniform(0, random.uniform(0,2))

    f = - (math.log(((1/0.995)-1),0.5))


    # initial condition 
    x0 = [0.1,1.7,0.1,0.05,0]

    
    # values of time 
    timespent = 10
    t = np.linspace(0,timespent,(timespent*50)+1) 

    p = np.linspace(0,timespent,(timespent*50)+1)
    n = np.linspace(0,timespent,(timespent*50)+1)
    u = np.linspace(0,timespent,(timespent*50)+1) 
    c = np.linspace(0,timespent,(timespent*50)+1) 
    s = np.linspace(0,timespent,(timespent*50)+1)

    p[0] = 0
    n[0] = 0.55
    u[0] = 1.7
    c[0] = 0.1
    s[0] = 0.05


    for i in range((timespent*50)):

        p[i+1] = p[i] + w * ((c[i]/(c[i]+1))) - p[i]*(u[i] * (s[i]/(ki+s[i]))*n[i]+dp)
        n[i+1] = n[i] + u[i] * (s[i]/(ki+s[i]))*n[i]
        u[i] = umax * (s[i]/(ks+s[i])) * ((c[i]/cmax)**(-f))/((0.5**(-f))+((c[i]/cmax)**(-f)))
        c[i+1] = c[i] + -c[i]*(u[i] * (s[i]/(ki+s[i]))*n[i]+dp)
        s[i+1] = s[i] + -(g) * u[i] * (s[i]/(ki+s[i]))*n[i]
    
    u[-1] = u[-2]

    # plt.plot(t,n, label="number of cells") 
    figure, axis = plt.subplots(2, 3) 
    

    axis[0, 0].plot(t, n) 
    axis[0, 0].set_title("number of cells") 

    axis[0, 1].plot(t, s) 
    axis[0, 1].set_title("Nutrition Concentration") 

    axis[0, 2].plot(t, u) 
    axis[0, 2].set_title("Growth Rate") 

    axis[1, 0].plot(t, p) 
    axis[1, 0].set_title("Protein Concentration") 

    axis[1, 1].plot(t, c) 
    axis[1, 1].set_title("Inducer Concentration") 


        
    plt.show()





