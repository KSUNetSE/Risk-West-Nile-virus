# -*- coding: utf-8 -*-
"""
Created on Thu Apr 18 12:45:47 2024

@author: shosseini
"""

import numpy as np
from scipy.integrate import odeint
from matplotlib import pyplot as plt
from scipy.optimize import curve_fit
from scipy.stats import norm
from scipy.stats import gamma
from scipy.stats import zscore
from scipy.stats import poisson
from scipy.stats import beta
import time

#--------------------------------------Section of importing data 
C2006 = list(np.loadtxt("./C2006.txt"))
T2006 = list(np.loadtxt("./T2006.txt"))
H2006 = list(np.loadtxt("./H2006.txt"))
#---------------------------------------------------------------
C2007 = list(np.loadtxt("./C2007.txt"))
T2007 = list(np.loadtxt("./T2007.txt"))
H2007 = list(np.loadtxt("./H2007.txt"))
#---------------------------------------------------------------
C2008 = list(np.loadtxt("./C2008.txt"))
T2008 = list(np.loadtxt("./T2008.txt"))
H2008 = list(np.loadtxt("./H2008.txt"))
#---------------------------------------------------------------
C2009 = list(np.loadtxt("./C2009.txt"))
T2009 = list(np.loadtxt("./T2009.txt"))
H2009 = list(np.loadtxt("./H2009.txt"))
#---------------------------------------------------------------
C2011 = list(np.loadtxt("./C2011.txt"))
T2011 = list(np.loadtxt("./T2011.txt"))
H2011 = list(np.loadtxt("./H2011.txt"))
#---------------------------------------------------------------
C2012 = list(np.loadtxt("./C2012.txt"))
T2012 = list(np.loadtxt("./T2012.txt"))
H2012 = list(np.loadtxt("./H2012.txt"))
#---------------------------------------------------------------
C2013 = list(np.loadtxt("./C2013.txt"))
T2013 = list(np.loadtxt("./T2013.txt"))
H2013 = list(np.loadtxt("./H2013.txt"))
#---------------------------------------------------------------
C2014 = list(np.loadtxt("./C2014.txt"))
T2014 = list(np.loadtxt("./T2014.txt"))
H2014 = list(np.loadtxt("./H2014.txt"))
#---------------------------------------------------------------
C2015 = list(np.loadtxt("./C2015.txt"))
T2015 = list(np.loadtxt("./T2015.txt"))
H2015 = list(np.loadtxt("./H2015.txt"))
#---------------------------------------------------------------
C2016 = list(np.loadtxt("./C2016.txt"))
T2016 = list(np.loadtxt("./T2016.txt"))
H2016 = list(np.loadtxt("./H2016.txt"))
#---------------------------------------------------------------
C2017 = list(np.loadtxt("./C2017.txt"))
T2017 = list(np.loadtxt("./T2017.txt"))
H2017 = list(np.loadtxt("./H2017.txt"))
#---------------------------------------------------------------
C2018 = list(np.loadtxt("./C2018.txt"))
T2018 = list(np.loadtxt("./T2018.txt"))
H2018 = list(np.loadtxt("./H2018.txt"))
#---------------------------------------------------------------
C2019 = list(np.loadtxt("./C2019.txt"))
T2019 = list(np.loadtxt("./T2019.txt"))
H2019 = list(np.loadtxt("./H2019.txt"))
#---------------------------------------------------------------
C2020 = list(np.loadtxt("./C2020.txt"))
T2020 = list(np.loadtxt("./T2020.txt"))
H2020 = list(np.loadtxt("./H2020.txt"))
#---------------------------------------------------------------
#C2021 = list(np.loadtxt("./C2021.txt"))
T2021 = list(np.loadtxt("./T2021.txt"))
H2021 = list(np.loadtxt("./H2021.txt"))

tt=0
Carying_2006=C2006[tt:365]
Temp_2006=T2006[tt:365]
Humidity_2006=H2006[tt:365]
#Precipitation_2006=P2006[tt:365]
#----------------------------------------------------------
Carying_2007=C2007[tt:365]
Temp_2007=T2007[tt:365]
Humidity_2007=H2007[tt:365]
#Precipitation_2007=P2007[tt:365]
#----------------------------------------------------------
Carying_2008=C2008[tt:365]
Temp_200=T2008[tt:365]
Humidity_2008=H2008[tt:365]
#Precipitation_2008=P2008[tt:365]
#----------------------------------------------------------
Carying_2009=C2009[tt:365]
Temp_2009=T2009[tt:365]
Humidity_2009=H2009[tt:365]
#Precipitation_2009=P2009[tt:365]
#----------------------------------------------------------
# Carying_2010=C2010[tt:365]
# Temp_2010=T2010[tt:365]
# Humidity_2010=H2010[tt:365]
#Precipitation_2010=P2010[tt:365]
#----------------------------------------------------------
Carying_2011=C2011[tt:365]
Temp_2011=T2011[tt:365]
Humidity_2011=H2011[tt:365]
#Precipitation_2011=P2011[tt:365]
#----------------------------------------------------------
Carying_2012=C2012[tt:365]
Temp_2012=T2012[tt:365]
Humidity_2012=H2012[tt:365]
#Precipitation_2012=P2012[tt:365]
#----------------------------------------------------------
Carying_2013=C2013[tt:365]
Temp_2013=T2013[tt:365]
Humidity_2013=H2013[tt:365]
#Precipitation_2013=P2013[tt:365]
#----------------------------------------------------------
Carying_2014=C2014[tt:365]
Temp_2014=T2014[tt:365]
Humidity_2014=H2014[tt:365]
#Precipitation_2014=P2014[tt:365]
#----------------------------------------------------------
Carying_2015=C2015[tt:365]
Temp_2015=T2015[tt:365]
Humidity_2015=H2015[tt:365]
#Precipitation_2015=P2015[tt:365]
#----------------------------------------------------------
Carying_2016=C2016[tt:365]
Temp_2016=T2016[tt:365]
Humidity_2016=H2016[tt:365]
#Precipitation_2016=P2016[tt:365]
#----------------------------------------------------------
Carying_2017=C2017[tt:365]
Temp_2017=T2017[tt:365]
Humidity_2017=H2017[tt:365]
#Precipitation_2017=P2017[tt:365]
#----------------------------------------------------------
Carying_2018=C2018[tt:365]
Temp_2018=T2018[tt:365]
Humidity_2018=H2018[tt:365]
#Precipitation_2018=P2018[tt:365]
#----------------------------------------------------------
Carying_2019=C2019[tt:365]
Temp_2019=T2019[tt:365]
Humidity_2019=H2019[tt:365]
#Precipitation_2019=P2019[tt:365]
#----------------------------------------------------------
Carying_2020=C2020[tt:365]
Temp_2020=T2020[tt:365]
Humidity_2020=H2020[tt:365]
#Precipitation_2020=P2020[tt:365]
#----------------------------------------------------------
# Carying_2021=C2021[tt:365]
# Temp_2021=T2021[tt:365]
# Humidity_2021=H2021[tt:365]
#Precipitation_2021=P2021[tt:365]
#----------------------------------------------------------
data_dict_Carying = {}
data_dict_Temp={}
data_dict_Humid={}
for h in range(0,25):
    diclist=[]
    diclisttem=[]
    diclistpre=[]
    for hum in Humidity_2006:
        if h==np.round(hum): 
            diclist.append(Carying_2006[(Humidity_2006.index(hum))])
            diclisttem.append(Temp_2006[(Humidity_2006.index(hum))])
            #diclistpre.append(Precipitation_2006[(Humidity_2006.index(hum))])
    for hum in Humidity_2007:
        if h==np.round(hum): 
            diclist.append(Carying_2007[(Humidity_2007.index(hum))])
            diclisttem.append(Temp_2007[(Humidity_2007.index(hum))])
            #diclistpre.append(Precipitation_2007[(Humidity_2007.index(hum))])
    # for hum in Humidity_2008:
    #     if h==np.round(hum): 
    #         diclist.append(Carying_2008[(Humidity_2008.index(hum))])
    #         diclisttem.append(Temp_2008[(Humidity_2008.index(hum))])
    #         #diclistpre.append(Precipitation_2008[(Humidity_2008.index(hum))])
    for hum in Humidity_2009:
        if h==np.round(hum): 
            diclist.append(Carying_2009[(Humidity_2009.index(hum))])
            diclisttem.append(Temp_2009[(Humidity_2009.index(hum))])
            #diclistpre.append(Precipitation_2009[(Humidity_2009.index(hum))])
    # for hum in Humidity_2010:
    #     if h==np.round(hum): 
    #         diclist.append(Carying_2010[(Humidity_2010.index(hum))])
    #         diclisttem.append(Temp_2010[(Humidity_2010.index(hum))])
    #         #diclistpre.append(Precipitation_2010[(Humidity_2010.index(hum))])
    for hum in Humidity_2011:
        if h==np.round(hum): 
            diclist.append(Carying_2011[(Humidity_2011.index(hum))])
            diclisttem.append(Temp_2011[(Humidity_2011.index(hum))])
            #diclistpre.append(Precipitation_2011[(Humidity_2011.index(hum))])
    for hum in Humidity_2012:
        if h==np.round(hum): 
            diclist.append(Carying_2012[(Humidity_2012.index(hum))])
            diclisttem.append(Temp_2012[(Humidity_2012.index(hum))])
            #diclistpre.append(Precipitation_2012[(Humidity_2012.index(hum))])
    for hum in Humidity_2013:
        if h==np.round(hum): 
            diclist.append(Carying_2013[(Humidity_2013.index(hum))])
            diclisttem.append(Temp_2013[(Humidity_2013.index(hum))])
            #diclistpre.append(Precipitation_2013[(Humidity_2013.index(hum))])
    for hum in Humidity_2014:
        if h==np.round(hum): 
            diclist.append(Carying_2014[(Humidity_2014.index(hum))])
            diclisttem.append(Temp_2014[(Humidity_2014.index(hum))])
            #diclistpre.append(Precipitation_2014[(Humidity_2014.index(hum))])
    for hum in Humidity_2015:
        if h==np.round(hum): 
            diclist.append(Carying_2015[(Humidity_2015.index(hum))])
            diclisttem.append(Temp_2015[(Humidity_2015.index(hum))])
            #diclistpre.append(Precipitation_2015[(Humidity_2015.index(hum))])
    for hum in Humidity_2016:
        if h==np.round(hum): 
            diclist.append(Carying_2016[(Humidity_2016.index(hum))])
            diclisttem.append(Temp_2016[(Humidity_2016.index(hum))])
            #diclistpre.append(Precipitation_2016[(Humidity_2016.index(hum))])
    for hum in Humidity_2017:
        if h==np.round(hum): 
            diclist.append(Carying_2017[(Humidity_2017.index(hum))])
            diclisttem.append(Temp_2017[(Humidity_2017.index(hum))])
            #diclistpre.append(Precipitation_2017[(Humidity_2017.index(hum))])
    for hum in Humidity_2018:
        if h==np.round(hum): 
            diclist.append(Carying_2018[(Humidity_2018.index(hum))])
            diclisttem.append(Temp_2018[(Humidity_2018.index(hum))])
            #diclistpre.append(Precipitation_2018[(Humidity_2018.index(hum))])
    for hum in Humidity_2019:
        if h==np.round(hum): 
            diclist.append(Carying_2019[(Humidity_2019.index(hum))])
            diclisttem.append(Temp_2019[(Humidity_2019.index(hum))])
            #diclistpre.append(Precipitation_2019[(Humidity_2019.index(hum))])
    for hum in Humidity_2020:
        if h==np.round(hum): 
            diclist.append(Carying_2020[(Humidity_2020.index(hum))])
            diclisttem.append(Temp_2020[(Humidity_2020.index(hum))])
            #diclistpre.append(Precipitation_2020[(Humidity_2020.index(hum))])
    if len(diclist)>0:
        data_dict_Carying[h]=list((diclist))
        data_dict_Temp[h]=list(diclisttem)


#---------------------------------------------------------------------------
def K_T_H(humidity,temperature):   
    rr=int(np.round(humidity))
    tem=temperature
    x = data_dict_Temp[rr]
    y = data_dict_Carying[rr]
    z_scores_x = np.abs((x - np.mean(x)) / np.std(x))
    z_scores_y = np.abs((y - np.mean(y)) / np.std(y))

# Define a threshold for outliers (you can adjust this threshold)
    outlier_threshold = 3

# Identify indices of outliers for both x and y
    outlier_indices = np.where((z_scores_y > outlier_threshold) | (z_scores_y > outlier_threshold))[0]
    x_no_outliers = np.delete(x, outlier_indices)
    y_no_outliers = np.delete(y, outlier_indices)
    poly_coefficients = np.polyfit(x_no_outliers, y_no_outliers, 10)
    x_range = np.linspace(min(x_no_outliers), max(x_no_outliers), 100)
    y_fit = np.abs(np.polyval(poly_coefficients, x_range))

    # print(min(x_no_outliers),max(x_no_outliers))
    if tem<min(x_no_outliers):
        car=np.abs(np.polyval(poly_coefficients, min(x_no_outliers)))
        # print("hhhhhhhhhhh",car)
    elif tem>min(x_no_outliers):
        car=np.abs(np.polyval(poly_coefficients, max(x_no_outliers)))
        # print("hhhhhhhhhhh",car)
    else:
        car=np.abs(np.polyval(poly_coefficients, tem))
        # print("hhhhhhhhhhh",car)

    # plt.scatter(x, y, label='Original Data')
    # plt.scatter(x_no_outliers, y_no_outliers, label='Data without Outliers')
    # plt.plot(x_range, y_fit, label='Fitted Polynomial (Order 11)', color='black')
    # plt.title("humidity range: {:.1f} - {:.1f}".format(rr-0.5, rr+0.5))
    # plt.xlabel('X-axis')
    # plt.ylabel('Y-axis')
    # plt.legend()
    # plt.show()
    return(car)
K_T_H(7,15)

#-----------------------------------function of carying capacity---------------
Hum=H2021[0:365]
dama=T2021[0:365]
Coematrix=np.zeros([9,365])
w=0

for t in range(0,365):
    if np.round(Hum[t])>max(data_dict_Carying.keys()):
        Coematrix[8,t]=K_T_H(max(data_dict_Carying.keys()),dama[t])-w
    elif np.round(Hum[t])<min(data_dict_Carying.keys()):
        Coematrix[8,t]=K_T_H(min(data_dict_Carying.keys()),dama[t])-w
    elif np.round(Hum[t]) not in data_dict_Carying.keys():
        number=np.round(Hum[t])
        my_vector=list(data_dict_Carying.keys())
        differences = np.abs(my_vector - number)
        closest_index = np.argmin(differences)
        closest_element = my_vector[closest_index]
        Coematrix[8,t]=K_T_H(closest_element,dama[t])-w
    else:
        Coematrix[8,t]=K_T_H(Hum[t],dama[t])-w
        
plt.plot(Coematrix[8,:],'-')
plt.show()


def tim2tem(time):
    tempc=dama[time]
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
 
for time in range(0,365):
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
    dacuMdt=Coematrix[5,time]*x[4]*max(0,(1-(x[5]/(Coematrix[8,time]))))-Coematrix[0,time]*Coematrix[1,time]*x[5]
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
t=np.linspace(0,365,364)
x0=(3000000,0,0,0,10,10,10,10,10,100,100,100,100,100,100)
x=odeint(BMPPT,x0,t)
R=x[0:365,3]
I=x[0:365,2]
M=x[0:365,6]
plt.plot(I,'--k',label='estimated profile of mosquitoes')

























