# -*- coding: utf-8 -*-
"""
Created on Thu Apr 25 12:58:31 2024

@author: shosseini
"""

import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from statsmodels.tsa.statespace.sarimax import SARIMAX
from statsmodels.tsa.arima.model import ARIMA
from statsmodels.tsa.ar_model import AutoReg
from scipy.integrate import odeint
from scipy.stats import gamma


# thid code will return the risk prediction for two weeks 
# the vector of the carying capacity (as well as temperature) should be importing to
# the code (so using another code the carying capacity should be predicted.)
# at the end it is plotted based in the perioi of iteration you want.

dama = list(np.loadtxt('2000-2021tem.txt'))
T=dama
prediction_lead=14
dama=dama[len(dama)-4*365:len(dama)]
for datenumber in range(100,365,1):
#for datenumber in {345}:
    # datenumber=260
    d=365-datenumber
    inuse=dama[0:len(dama)-d]
    #print(inuse)
    model = AutoReg(inuse, lags=365)
    model_fitted = model.fit()
    predictions = model_fitted.predict(start=len(inuse), end=len(inuse)+prediction_lead)
    T2=list(predictions)
    
    # plt.plot([i for i in range(365-d,365-d+16)],predictions,'b-',label='Prediction lead time 15 days')
    # plt.plot(dama[len(dama)-365:len(dama)-d],':r')
    # plt.xlabel('day')
    # plt.ylabel('temp')
    # plt.title('Short-Term Prediction of Temperature in Orange County, 2021')
    # plt.legend()
    # plt.savefig(f'forecast-short{datenumber}'+'.pdf')
    # plt.show()
    
    temp=T[len(T)-365:len(T)-365+datenumber]
    for i in range(0,len(predictions)):
        temp.append(np.round(predictions[i],decimals=2))
    # print(len(temp))
    
    #------------------------------------part two
    def tim2tem(time):
        tempc=temp[time]
        #if tempc>34:tempc=34
        # tim2tem=(t2t-32)*(5/9)
        return(tempc)
    def gammaMP(time):#mosquitoe development rate*****Briere********MDR*************
        t=time
        if tim2tem(t)<0.1 or tim2tem(t)>38.5:gammma=0
        elif tim2tem(t)>=0.1 and tim2tem(t)<=38.5:
            gammma=tim2tem(t)*(3.76*10**(-5))*(tim2tem(t)-0.1)*(np.sqrt(38.5-tim2tem(t)))
        return(gammma)
    def pLSP(time):#proportion of larval survival*****quadratic*****pLA*************
        t=time
        if tim2tem(t)<7.8 or tim2tem(t)>38.4:
            PLS=0
        elif tim2tem(t)>=7.8 and tim2tem(t)<38.4:
            PLS=(-3.6*(10**(-3)))*(tim2tem(t)-7.8)*(tim2tem(t)-38.4)
        return(PLS)
    def muMP(time): #mortality rate*******************linear********muM*************
        t=time
        if tim2tem(t)>=34.93 :mu=80000
        elif(tim2tem(t)<14):mu=6/99
        elif tim2tem(t)>=14 and tim2tem(t)<34.93:
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
    
    # we will fill this matrix by coeficient of each day so we have 365 set of data
    #------------------------------------------------------------------------------
    Coematrix=np.zeros([9,datenumber+2+prediction_lead])    
    N=[]
    
    
    
    #-----------------------------ploting profile of the mosquitoe for each year---
    #-------------Also it collects the sample of number of effective populations---
    #------------------------------------------------------------------------------
    
    
    C = list(np.loadtxt("./C2021.txt"))
    Casetime=list(np.loadtxt(f"./orangecases2021.txt"))
    def nonzero_indices(lst):
        return [index for index, value in enumerate(lst) if value != 0]
    ct=nonzero_indices(Casetime)
    #print(ct)
    for time in range(0,datenumber+1+prediction_lead):
        Coematrix[0,time]=gammaMP(time)
        Coematrix[1,time]=pLSP(time)
        Coematrix[2,time]=muMP(time)
        Coematrix[3,time]=BRP(time)
        Coematrix[4,time]=alphaEP(time)
        Coematrix[5,time]=bEP(time)
        Coematrix[6,time]=PDRP(time)
        Coematrix[7,time]=VCP(time)
    
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
        muHuman=(1/(77*365))
        #***** equations related to dynamic of Culex pipiens***********************
        #kC=D
        #--------------------Human equations
        CB=20000
        TH=x[0]+x[1]+x[2]+x[3]
        D=x[11]+x[12]+x[13]+x[14]
        ph=Coematrix[3,time]/(D+TH)
        dsHdt=(muHuman*(TH))-ph*Coematrix[7,time]*x[8]*x[0]-muHuman*x[0]
        #dsHdt=35-ph*Coematrix[7,time]*x[8]*x[0]-muHuman*x[0]
        deHdt=ph*Coematrix[7,time]*x[8]*x[0]-(incHumanWND)*x[1]-muHuman*x[1]
        diHdt=(incHumanWND)*x[1]-(infectionhumanWND)*x[2]-(muhumanWND)*x[2]-muHuman*x[2]
        drHdt=(infectionhumanWND)*x[2]-muHuman*x[3]
        #--------------------Mosquitoes' developement equations
        deggMdt=Coematrix[4,time]*(x[6]+x[7]+x[8])-Coematrix[5,time]*x[4]
        dacuMdt=Coematrix[5,time]*x[4]*max(0,(1-(x[5]/(C[time]))))-Coematrix[0,time]*Coematrix[1,time]*x[5]
        #-------------------Mosquitoe +pathogen+Bird+Human
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
    t=np.linspace(0,datenumber+1+prediction_lead,datenumber+prediction_lead)
    x0=(3000000,0,0,0,10,10,10,10,10,100,100,100,100,100,100)
    x=odeint(BMPPT,x0,t)
    R=x[0:datenumber+1+prediction_lead,3]
    I=x[0:datenumber+1+prediction_lead,2]
    M=x[0:datenumber+1+prediction_lead,6]
    #plt.plot(M,'--k',label='estimated profile of mosquitoes')
    # for i in ct:
    #     if i==min(ct):
    #         plt.axvline(x=i*7, color='b', linestyle='--', label='spillover time')
    #     else:
    #         plt.axvline(x=i*7, color='b', linestyle='--')
    
    # plt.xlabel('days')
    # plt.ylabel('Estimated number of mosquitoes')
    # plt.title('Orange County - prediction for 15 days, 2021, day'+str(datenumber))
    # plt.legend()
    
    # plt.savefig(f'Profile_orange_county{year}.jpg', format='jpg', dpi=300)
    # plt.savefig(f'Profile_orange_county{year}'+'.pdf')
    # plt.show()
    a=6.5888138464387
    loc=0
    scale=3429940769.08221
    probabilities = [0.03, 0.3, 0.6]
    quantiles = gamma.ppf(probabilities, a, loc, scale)
    
    for i in range(0,datenumber+prediction_lead):
        if M[i]<quantiles[0]:
            plt.axvline(i, color='g')
        elif M[i]>quantiles[0] and M[i]<quantiles[1]:
            plt.axvline(i, color='y')
        elif M[i]>quantiles[1] and M[i]<quantiles[2]:
            plt.axvline(i, color='orange')
        elif M[i]>quantiles[2]:
            plt.axvline(i, color='red')
    # plt.title('Probability of spillover for 200'+str(2021))
    # plt.xlabel('day')
    # plt.legend(['green:P < 0.03', 'yellow: 0.03 < P < 0.3', 'orange: 0.3 < P < 0.6', 'red: 0.6 < P '],
    #           loc='upper left')
    #plt.savefig('200' + str(f'{r}') + '.pdf')
    # plt.savefig('200' + str(f'{year}'))
    plt.ylim(0, 11.5**10)
    plt.xlim(0,370)
    plt.yticks([])
    plt.plot(M,'--w',label='Estimated profile of mosquitoes')
    #plt.axvline(s, color='b', linestyle='--', label='spillover time')

    #print(ct*7)
    for i in ct:
        if i==min(ct):
            plt.axvline(x=i*7, color='b', linestyle='--', label='spillover time')
        else:
            plt.axvline(x=i*7, color='b', linestyle='--')
    #plt.savefig(f'risk{datenumber}.jpeg', format='JPEG', dpi=500)
    

    plt.title('Orange County, Probabilistic Risk of Spillover for 20'+str(21))
    plt.xlabel('Days of the year')
    #plt.ylabel('Mosquito Profile')
    # plt.legend(['green:P < 0.03', 'yellow: 0.03 < P < 0.3', 'orange: 0.3 < P < 0.6', 'red: 0.6 < P '],
    #           loc='upper left')
    plt.axvline(x=datenumber, color='k', linestyle='--', linewidth=1, label='prediction period')
    plt.axvline(x=datenumber+prediction_lead, color='k', linestyle='--', linewidth=1)
    #plt.legend()
    plt.legend(loc='upper left', fontsize='small')  # You can also use specific sizes like 10, 12, etc.
    # if datenumber in {345}:
    plt.savefig(f'pic{datenumber}.pdf', format='pdf')

    plt.show()
    
    
