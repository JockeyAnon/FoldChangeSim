
import numpy as np
import math 
from scipy.integrate import odeint 
import random
import matplotlib.pyplot as plt 


newdir = "C:/Users/Daniel Park/Desktop/SimGraphs/"

for i in range(300):
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

    def odes(x,t):

        n=x[0]
        u=x[1]
        c=x[2]
        s=x[3]
        p=x[4]

        dndt = u * (s/(ki+s))*n

        dudt = umax * (s/(ks+s)) * ((c/cmax)**(-f))/((0.5**(-f))+((c/cmax)**(-f)))

        # dcdt = -c*(u+(qmax/(c+ks+((c**2)/ki))))
        dcdt = -c*(u+dp)

        dsdt = -(g) * dndt

        # dpdt = w * ((c/(c+1))) - p*(u+(qmax/(p+ks+((p**2)/ki))))
        dpdt = w * ((c/(c+1))) - p*(u+dp)

        return [dndt,dudt,dcdt,dsdt,dpdt]

    # initial condition 
    x0 = [0.1,1.7,0.1,0.05,0]

    
    # values of time 
    timespent = 10
    t = np.linspace(0,timespent,(timespent*50)+1) 
    
    # solving ODE 
    x = odeint(odes,x0,t) 
    n=x[:,0]
    u=x[:,1]
    c=x[:,2]
    s=x[:,3]
    p=x[:,4]

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


        
    filename = newdir + 'sim pic with variables '+ "ki = "+ str(ki) +"ks = "+ str(ks) + "dc = "+ str(dc) + "g = "+ str(g) +"dp = "+ str(dp) + "qmax = "+ str(qmax) + '.png'
    figure.savefig(filename)
    plt.close(figure)





