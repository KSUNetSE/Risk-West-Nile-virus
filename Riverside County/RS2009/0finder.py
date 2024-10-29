# -*- coding: utf-8 -*-
"""
Created on Fri Apr  5 11:30:14 2024

@author: shosseini
"""

import numpy as np
from scipy.optimize import curve_fit
import matplotlib.pyplot as plt
from scipy.integrate import odeint
from numba import jit
from scipy.integrate import solve_ivp
def logistic_function(t, L, k, t0):
    return L / (1 + np.exp(-k * (t - t0)))
t = np.linspace(0, 52, 52)  # Time vector

values=np.cumsum(list(np.loadtxt("./RS2009.txt")))
dama = np.loadtxt('RS-tem-2009.txt')
# plt.plot(dama)
# plt.show()




popt, pcov = curve_fit(logistic_function, t, values, p0=[max(values), 1, np.median(t)])
print(f"Optimized parameters: L={popt[0]:.2f}, k={popt[1]:.2f}, t0={popt[2]:.2f}")

tt = np.linspace(0, 52, 364)
t=[i for i in range(0,364)]

y=logistic_function(tt, *popt)




exp_data=y

def tim2tem(time):
    tempc=dama[time]
    #if tempc>34.93:tempc=34.9
    return(tempc)
def gammaMP(time):#mosquitoe development rate*****Briere********MDR*************
    t=time
    if tim2tem(t)<0.1 or tim2tem(t)>38.5:gamma=0
    elif tim2tem(t)>=0.1 and tim2tem(t)<=38.5:
        gamma=tim2tem(t)*(3.76*10**(-5))*(tim2tem(t)-0.1)*(np.sqrt(38.5-tim2tem(t)))
    return(gamma)
def pLSP(time):#proportion of larval survival*****quadratic*****pLA*************
    t=time
    if tim2tem(t)<7.8 or tim2tem(t)>38.4:
        PLS=0
    elif tim2tem(t)>=7.8 and tim2tem(t)<38.4:
        PLS=(-3.6*(10**(-3)))*(tim2tem(t)-7.8)*(tim2tem(t)-38.4)
    return(PLS)
def muMP(time): #mortality rate*******************linear********muM*************
    t=time
    if tim2tem(t)>34.93 :
        mu=123456789
    if(tim2tem(t)<14):
        mu=6/99
    elif tim2tem(t)>=14 and tim2tem(t)<=34.93:
        mu=(-4.86*tim2tem(t)+169.8)**(-1)
    return(mu)
def BRP(time):# a***biting rate*******************Briere********BR**************
    t=time
    if tim2tem(t)<9.4 or tim2tem(t)>39.6:
        Br=0
    elif tim2tem(t)>=9.4 and tim2tem(t)<39.6:
        Br=tim2tem(t)*(1.7*10**(-4))*(tim2tem(t)-9.4)*(np.sqrt(39.6-tim2tem(t)))
    return(Br)
def alphaEP(time):#fecundity**********************quadratic*****EFGC/LGC********
    t=time
    if tim2tem(t)<5.3 or tim2tem(t)>38.9:alphaE=0
    elif tim2tem(t)>=5.3 and tim2tem(t)<38.9:
        alphaE=(-(5.98*10**(-1))*(tim2tem(t)-5.3)*(tim2tem(t)-38.9))*BRP(t)
    return(alphaE)
def bEP(time):#egg viability**********************quadratic*****EV**************
    t=time
    if tim2tem(t)<3.2 or tim2tem(t)>42.6:bE=0
    elif tim2tem(t)>=3.2 and tim2tem(t)<42.6:
        bE=-(2.11*10**(-3))*(tim2tem(t)-3.2)*(tim2tem(t)-42.6)
    return(bE)
def PDRP(time):#*********briere function***pathogen Development Rate************
    t=time
    if tim2tem(t)<=11.4 or tim2tem(t)>45.2:pdr=0
    elif tim2tem(t)>11.4 and tim2tem(t)<=45.2:
        pdr=tim2tem(t)*(7.38*10**(-5))*(tim2tem(t)-11.4)*(np.sqrt(45.2-tim2tem(t)))
    return(pdr)
def VCP(time):#instead od trasmission efficiency
    t=time
    if tim2tem(t)<16.8 or tim2tem(t)>38.9:v=0
    elif tim2tem(t)>=16.8 and tim2tem(t)<=38.9:
        v=-(3.05*10**(-3))*(tim2tem(t)-16.8)*(tim2tem(t)-38.9)
    return(v)
Coematrix=np.zeros([11,365])    #we will fill this matrix by coeficient of each day so we have 365 set of data
for time in range(0,365):
    Coematrix[0,time]=gammaMP(time)
    Coematrix[1,time]=pLSP(time)
    Coematrix[2,time]=muMP(time)
    Coematrix[3,time]=BRP(time)
    Coematrix[4,time]=alphaEP(time)
    Coematrix[5,time]=bEP(time)
    Coematrix[6,time]=PDRP(time)
    Coematrix[7,time]=VCP(time)


C=516000000
@jit
def BMPPT(x,time):
    time=int(divmod(time, 365)[1])
    x[0]=max(0,x[0]);x[1]=max(0,x[1]);x[2]=max(0,x[2]);x[3]=max(0,x[3]);x[4]=max(0,x[4])
    x[5]=max(0,x[5]);x[6]=max(0,x[6]);x[7]=max(0,x[7]);x[8]=max(0,x[8]);x[9]=max(0,x[9]);
    x[11]=max(0,x[11]);x[12]=max(0,x[12]);x[13]=max(0,x[13]);x[14]=max(0,x[14])
    #**************************************************************************
    sigma=15;mu=110
    #FB=1.96*(1/np.sqrt((2*np.pi)*(sigma**2)))*np.exp(-(2*sigma**(-2))*(((np.divmod(t,365)[1])-mu)**2))
    FB=2
    EVB=0.021
    muE=4.9315*(10**(-4))
    BDR=9.041*(10**(-4))
    muF=6.30136*(10**(-4))
    muB=(0.15/365)
    incB=1/3
    recoveryrateC=1/6
    muwnvB=0.9
    muHbirth=5
    muhumanWND=0.01
    incHumanWND=1/6
    infectionhumanWND=1/3
    #muHuman=(1/(77*365))
    muHuman=0.00002137
    #***** equations related to dynamic of Culex pipiens***********************
    CB=1217596
    #--------------------Human equations
    TH=x[0]+x[1]+x[2]+x[3]
    D=x[11]+x[12]+x[13]+x[14]
    ph=Coematrix[3,time]/(D+TH)
    dsHdt=muHuman-ph*Coematrix[7,time]*x[8]*x[0]-muHuman*x[0]
    #dsHdt=35-ph*Coematrix[7,time]*x[8]*x[0]-muHuman*x[0]
    deHdt=ph*Coematrix[7,time]*x[8]*x[0]-(incHumanWND)*x[1]-muHuman*x[1]
    diHdt=(incHumanWND)*x[1]-(infectionhumanWND)*x[2]-(muhumanWND)*x[2]-muHuman*x[2]
    drHdt=(infectionhumanWND)*x[2]-muHuman*x[3]
    #--------------------Mosquitoes' developement equations
    deggMdt=Coematrix[4,time]*(x[6]+x[7]+x[8])-Coematrix[5,time]*x[4]

    dacuMdt=Coematrix[5,time]*x[4]*max(0,(1-(x[5]/C)))-Coematrix[0,time]*Coematrix[1,time]*x[5]
    #-------------------Mosquitoe +pathogen+Bird+Human
    # dsMdt=Coematrix[0,time]*Coematrix[1,time]*x[5]-ph*x[13]*x[6]-Coematrix[2,time]*x[6]
    # deMdt=ph*x[13]*x[6]-Coematrix[6,time]*x[7]-Coematrix[2,time]*x[7]
    # diMdt=Coematrix[6,time]*x[7]-Coematrix[2,time]*x[8]
    if Coematrix[2,time]==123456789:
        dsMdt=0
        deMdt=0
        diMdt=0
    else:
        dsMdt=Coematrix[0,time]*Coematrix[1,time]*x[5]-ph*x[13]*x[6]-Coematrix[2,time]*x[6]
        deMdt=ph*x[13]*x[6]-Coematrix[6,time]*x[7]-Coematrix[2,time]*x[7]
        diMdt=Coematrix[6,time]*x[7]-Coematrix[2,time]*x[8]
    #------------------Birds' developement equations
    deggBdt=FB*(x[11]+x[12]+x[13]+x[14])-EVB*x[9]-muE*x[9]
    dfleBdt=EVB*x[9]*max(0,(1-(x[10]/CB)))-BDR*x[10]-muF*x[10]
    #------------------Bird+pathogen+Mosquito
    dsBdt=BDR*x[10]-ph*Coematrix[7,time]*x[8]*x[11]-muB*x[11]
    deBdt=ph*Coematrix[7,time]*x[8]*x[11]-incB*x[12]-muB*x[12]
    diBdt=incB*x[12]-muwnvB*x[13]-muB*x[13]-recoveryrateC*x[13]
    drBdt=recoveryrateC*x[13]-muB*x[14]
    dxdt=(dsHdt,deHdt,diHdt,drHdt,deggMdt,dacuMdt,dsMdt,deMdt,diMdt,deggBdt,dfleBdt,dsBdt,deBdt,diBdt,drBdt)
    return(dxdt)
r=2
d=r-1
dd=r
t=np.linspace(0,r*365,r*364)
x0=(2474000,0,0,0,1,1,1,1,0,100,100,100,100,100,100)
x=odeint(BMPPT,x0,t)
MI=x[0:r*365,8]
MS=x[0:r*365,6]
R=x[0:r*365,3]

fig, axs = plt.subplots(3, 1, figsize=(10, 10), sharex=True)

# Plot R
axs[0].plot(MI[(d)*365:dd*365], label='MI')
axs[0].set_ylabel('R')
axs[0].legend()

# Plot I
axs[1].plot(R[(d)*365:dd*365], label='Removed', color='orange')
axs[1].plot(logistic_function(tt, *popt), label='real cases', color='red')
axs[1].set_ylabel('I')
axs[1].set_xlabel('Time')
axs[1].legend()

# Plot 
axs[2].plot(MS[(d)*365:dd*365], label='MS', color='orange')
axs[2].set_ylabel('I')
axs[2].set_xlabel('Time')
axs[2].legend()
plt.tight_layout()
plt.show()
