# -*- coding: utf-8 -*-
"""
Created on Tue Apr 16 11:26:55 2024

@author: shosseini
"""

import numpy as np
from scipy.integrate import odeint
from matplotlib import pyplot as plt
from scipy.stats import gamma
from scipy.optimize import minimize
import pandas as pd
from scipy.stats import norm, multivariate_normal, chi2
from numba import jit
#----------------------------------------------------------------
# This code retrieves information such as temperature, number of cases, and 
# carrying capacity for various years. 
# It performs two main functions:
# 1. Simultaneously plots the mosquito profile (or any other required data)
# along with the times of spillover events.
# 2. Provides an estimation of risk for each year based on defined quintiles
# (we have set these at 3%, 30%, and 60%).
# Additionally, it allows for simultaneous examination of spillover timing 
# and probabilistic risk assessment. Plus ploting the pdf of spillover

def get_temperature_for_year(file_path, year):
    # Read the CSV file
    df = pd.read_csv(file_path)
    df_year = df[df['Year'] == year]
    # Convert the temperature column to a list of real values (floats)
    t = df_year['T'].tolist()
    h = df_year['H'].tolist()
    p = df_year['P'].tolist()
    return t,h,p
file_path = 'Riversideinfo.csv'



T2000 =get_temperature_for_year(file_path, 2000)[0]
T2001 =get_temperature_for_year(file_path, 2001)[0]
T2002 =get_temperature_for_year(file_path, 2002)[0]
T2003 =get_temperature_for_year(file_path,2003)[0]
T2004 =get_temperature_for_year(file_path, 2004)[0]
T2005 = get_temperature_for_year(file_path, 2005)[0]
#-------------------------------------------------------------------------
T2006 = get_temperature_for_year(file_path, 2006)[0]
T2007 = get_temperature_for_year(file_path, 2007)[0]
T2008 = get_temperature_for_year(file_path, 2008)[0]
T2009 = get_temperature_for_year(file_path, 2009)[0]
T2010 = get_temperature_for_year(file_path, 2010)[0]
T2011 = get_temperature_for_year(file_path, 2011)[0]
T2012 = get_temperature_for_year(file_path, 2012)[0]
T2013 = get_temperature_for_year(file_path, 2013)[0]
T2014 = get_temperature_for_year(file_path, 2014)[0]
T2015 = get_temperature_for_year(file_path, 2015)[0]
T2016 = get_temperature_for_year(file_path, 2016)[0]
T2017 = get_temperature_for_year(file_path, 2017)[0]
T2018 = get_temperature_for_year(file_path, 2018)[0]
T2019 = get_temperature_for_year(file_path, 2019)[0]
T2020 = get_temperature_for_year(file_path, 2020)[0]
#-------------------------------------------------------------------------
# T2021 = get_temperature_for_year(file_path, 2021)[0]
# T2022 = get_temperature_for_year(file_path, 2022)[0]
# T2023 = get_temperature_for_year(file_path, 2023)[0]



def moving_average_numpy_full(data, window_size):
    """Compute the moving average using NumPy with padding."""
    window = np.ones(int(window_size)) / float(window_size)
    extended_data = np.pad(data, pad_width=(window_size//2, window_size-1-window_size//2), mode='edge')
    return np.convolve(extended_data, window, 'valid')
window_size=1


#--------------------------------------Section of importing data 
C2006 = list(np.loadtxt("./C2006.txt"))
C2007 = list(np.loadtxt("./C2007.txt"))
C2008 = list(np.loadtxt("./C2008.txt"))
C2009 = list(np.loadtxt("./C2009.txt"))
C2010=C2009
C2011 = list(np.loadtxt("./C2011.txt"))
C2012 = list(np.loadtxt("./C2012.txt"))
C2013 = list(np.loadtxt("./C2013.txt"))
C2014 = list(np.loadtxt("./C2014.txt"))
C2015 = list(np.loadtxt("./C2015.txt"))
C2016 = list(np.loadtxt("./C2016.txt"))
C2017 = list(np.loadtxt("./C2017.txt"))
C2018 = list(np.loadtxt("./C2018.txt"))
C2019 = list(np.loadtxt("./C2019.txt"))
C2020 = list(np.loadtxt("./C2020.txt"))
C2006 = moving_average_numpy_full(C2006, window_size)
C2007 = moving_average_numpy_full(C2007, window_size)
C2008 = moving_average_numpy_full(C2008, window_size)
C2009 = moving_average_numpy_full(C2009, window_size)
C2010 = C2009
C2011 = moving_average_numpy_full(C2011, window_size)
C2012 = moving_average_numpy_full(C2012, window_size)
C2013 = moving_average_numpy_full(C2013, window_size)
C2014 = moving_average_numpy_full(C2014, window_size)
C2015 = moving_average_numpy_full(C2015, window_size)
C2016 = moving_average_numpy_full(C2016, window_size)
C2017 = moving_average_numpy_full(C2017, window_size)
C2018 = moving_average_numpy_full(C2018, window_size)
C2019 = moving_average_numpy_full(C2019, window_size)
C2020 = moving_average_numpy_full(C2020, window_size)

#---------------------------------------------------------------
T2021 = list(np.loadtxt("./2021forecasted_tem.txt"))
C2021 = list(np.loadtxt("./C_esti_2021.txt"))
#---------------------------------------------------------------
#---------------------------------------------------------------
T2022 = list(np.loadtxt("./2022forecasted_tem.txt"))
C2022 = list(np.loadtxt("./C_esti_2022.txt"))
#---------------------------------------------------------------
T2023 = list(np.loadtxt("./2023forecasted_tem.txt"))
C2023 = list(np.loadtxt("./C_esti_2023.txt"))
#--------------------------------------Section of the coeficient of
#                               the nonlinear model for pipiens mosquitoes----------------


def tim2tem(time):
    tempc=dama[time]
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
    elif tim2tem(t)>=5.3 and tim2tem(t)<=38.9:
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
Coematrix=np.zeros([9,365])    
N=[]
BR=[]
#-----------------------------ploting profile of the mosquitoe for each year---
#-------------Also it collects the sample of number of effective populations---
#------------------------------------------------------------------------------

#years = [year for year in range(2006, 2020) if year !=2010]
years = [year for year in range(2006, 2021) if year !=2010]
for year in years:
    redyear=[]
    oranyear=[]
    yellyear=[]
    redtogreen=[]
    greenyear=[]
    risktoyear=[]
    relativeredtorisk=[]
    dama = list(get_temperature_for_year(file_path, year)[0])
    if year == 2010:
        C=list(np.loadtxt("./C2009.txt"))
    else:
        C=list(np.loadtxt(f"./C{year}.txt"))
    Casetime=list(np.loadtxt(f"./RS{year}.txt"))
    def nonzero_indices(lst):
        return [index for index, value in enumerate(lst) if value != 0]
    ct=nonzero_indices(Casetime)
    #print(ct)
    for time in range(0,365):
        Coematrix[0,time]=gammaMP(time)
        Coematrix[1,time]=pLSP(time)
        Coematrix[2,time]=muMP(time)
        Coematrix[3,time]=BRP(time)
        Coematrix[4,time]=alphaEP(time)
        Coematrix[5,time]=bEP(time)
        Coematrix[6,time]=PDRP(time)
        Coematrix[7,time]=VCP(time)
    #@jit
    def BMPPT(x,time):
        time=int(divmod(time, 365)[1])
        x[0]=max(0,x[0]);x[1]=max(0,x[1]);x[2]=max(0,x[2]);x[3]=max(0,x[3]);x[4]=max(0,x[4])
        x[5]=max(0,x[5]);x[6]=max(0,x[6]);x[7]=max(0,x[7]);x[8]=max(0,x[8]);x[9]=max(0,x[9]);
        x[11]=max(0,x[11]);x[12]=max(0,x[12]);x[13]=max(0,x[13]);x[14]=max(0,x[14])
        #**************************************************************************
        FB=2
        EVB=0.021
        muE=4.9315*(10**(-4))
        BDR=9.041*(10**(-4))
        muF=6.30136*(10**(-4))
        muB=(0.15/365)
        incB=1/3
        recoveryrateC=1/6
        muwnvB=0.9
    
        muhumanWND=0.01
        incHumanWND=1/6
        infectionhumanWND=1/3
        #muHuman=(1/(77*365))
        muHuman=0.00002137
        #***** equations related to dynamic of Culex pipiens***********************
        #--------------------Human equations
        CB=1217596
        TH=x[0]+x[1]+x[2]+x[3]
        D=x[11]+x[12]+x[13]+x[14]
        ph=Coematrix[3,time]/(D+TH)
        dsHdt=(muHuman)-ph*Coematrix[7,time]*x[8]*x[0]-muHuman*x[0]
        deHdt=ph*Coematrix[7,time]*x[8]*x[0]-(incHumanWND)*x[1]-muHuman*x[1]
        diHdt=(incHumanWND)*x[1]-(infectionhumanWND)*x[2]-(muhumanWND)*x[2]-muHuman*x[2]
        drHdt=(infectionhumanWND)*x[2]-muHuman*x[3]
        #--------------------Mosquitoes' developement equations
        deggMdt=Coematrix[4,time]*(x[6]+x[7]+x[8])-Coematrix[5,time]*x[4]
    
        dacuMdt=Coematrix[5,time]*x[4]*max(0,(1-(x[5]/(C[time]))))-Coematrix[0,time]*Coematrix[1,time]*x[5]
        #-------------------Mosquitoe +pathogen+Bird+Human
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
    x0=(2474000,0,0,0,1,1,1,1,0,100,100,100,100,100,100)
    r=2
    t=np.linspace(0,r*365,r*364)
    x=odeint(BMPPT,x0,t)
    R=x[(r-1)*365-1:r*365,3]
    I=x[(r-1)*365-1:r*365,2]
    M=x[(r-1)*365-1:r*365,6]
    
    plt.plot(M/(max(M)),'--k',label='Estimated profile of mosquitoes')
    #plt.plot(Coematrix[3,:]/max(Coematrix[3,:]),'b:')
    # plt.plot(M,'--k',label='Estimated profile of mosquitoes')
    for i in ct:
        if i==min(ct):
            plt.axvline(x=i*7, color='r', linestyle='--', label='Spillover time')
        else:
            plt.axvline(x=i*7, color='r', linestyle='--')
    plt.xlabel('days')
    plt.ylabel('Estimated Relative Abundance of Mosquitoes')
    plt.title('Riverside County - '+str(year))
    plt.legend()
    plt.savefig(f'Profile_riverside_county{year}.jpg', format='jpg', dpi=1000)
    # plt.savefig(f'Profile_Riverside_county{year}'+'.pdf')
    # plt.figure(figsize=(10, 6))
    plt.legend()
    plt.show()
    BR.append(Coematrix[3,ct[0]*7])
    N.append(M[ct[0]*7])
    # for i in ct:
    #     if Coematrix[3,i*7]!=0:
    #         BR.append(Coematrix[3,i*7])
    #         N.append(M[i*7])
        
print("BR",BR)
print("N",N)




BR_array = np.array(BR)
N_array = np.array(N)
def fit_gamma(data):
    def negative_log_likelihood(params):
        shape, scale = params
        return -np.sum(gamma.logpdf(data, shape, scale=scale))
    
    initial_params = [1.0, np.mean(data)]
    result = minimize(negative_log_likelihood, initial_params, bounds=[(0, None), (0, None)])
    return result.x
params_N = fit_gamma(N_array)
shape_N, scale_N = params_N

params_BR = fit_gamma(BR_array)
shape_BR, scale_BR = params_BR

F_N = gamma.cdf(N_array, shape_N, scale=scale_N)
F_BR = gamma.cdf(BR_array, shape_BR, scale=scale_BR)

U = norm.ppf(F_N)
V = norm.ppf(F_BR)

correlation_coefficient = np.corrcoef(BR_array, N_array)[0, 1]
mean = [0, 0]
cov = [[1, correlation_coefficient], [correlation_coefficient, 1]]
rv = multivariate_normal(mean, cov)

n = np.linspace(0, 5*max(N_array), 100)
br = np.linspace(0, 5*max(BR_array), 100)
N_grid, BR_grid = np.meshgrid(n, br)

U_grid = norm.ppf(gamma.cdf(N_grid, shape_N, scale=scale_N))
V_grid = norm.ppf(gamma.cdf(BR_grid, shape_BR, scale=scale_BR))
pos = np.dstack((U_grid, V_grid))
Z = rv.pdf(pos)
probability_levels = [0.6, 0.9, 0.95]
chi2_values = [chi2.ppf(level, df=2) for level in probability_levels]
contour_levels = [rv.pdf([0, 0]) * np.exp(-0.5 * chi2_val) for chi2_val in chi2_values]
fig, ax = plt.subplots(figsize=(12, 8))
contour_filled = ax.contourf(N_grid, BR_grid, Z, levels=50, cmap='viridis')
plt.colorbar(contour_filled, ax=ax, label='Joint PDF')
contour_lines = plt.contour(N_grid, BR_grid, Z, levels=sorted(contour_levels), colors=['red', 'orange', 'black'])
ax.clabel(contour_lines,  fmt={contour_levels[0]: probability_levels[0] , contour_levels[1]: probability_levels[1], contour_levels[2]: probability_levels[1]}, inline=True, fontsize=10)
ax.set_xlabel('Abundance of Mosquitoes')
ax.set_ylabel('Biting Rate')
ax.set_title('Joint PDF with Contour Levels')
plt.grid(True)
plt.savefig(f'Joint PDF with Contour Levels{year}-Riverside.jpg', format='jpg', dpi=1000)
plt.show()


# x0 = 600000000  # Example value for N
# y0 = 0.24 
# plt.plot(x0, y0, 'r*')
# plt.show()
# i=0

# A = contour_lines.collections[0].get_paths()
# for path in A:  # Loop through each path
#     if path.contains_point((x0, y0)):  # Check if the point is within the path
#         i=i+1
#         break  # If the point is found in any path, you may want to break the loop
# A = contour_lines.collections[1].get_paths()
# for path in A:  # Loop through each path
#     if path.contains_point((x0, y0)):  # Check if the point is within the path
#         i=i+1
#         break  # If the point is found in any path, you may want to break the loop
# A = contour_lines.collections[2].get_paths()
# for path in A:  # Loop through each path
#     if path.contains_point((x0, y0)):  # Check if the point is within the path
#         i=i+1
#         break  # If the point is found in any path, you may want to break the loop
# print(i)
# if i==1:print("yellow") 
# elif i==2:print("orange")
# elif i==3: print("red")
# elif i==0:print("green") 







#------------------------------------------------------------------------------
#------------------section of plotting the pdf of spillover--------------------
#------------------------------------------------------------------------------

a, loc, scale = gamma.fit(N, floc=0)  # We fix the location to zero
# print(a, loc, scale)
x = np.linspace(0, max(N), 1000)
pdf_fitted = gamma.pdf(x, a, loc, scale)
plt.figure(figsize=(10, 6))
plt.hist(N, bins=5, density=True, color='blue', edgecolor='black', alpha=0.6)
plt.plot(x, pdf_fitted, 'r-', label=f'Gamma PDF\nshape={a:.2f}, scale={scale:.2f}')
plt.title('Histogram of Data with Fitted Gamma PDF, Riverside County, California')
plt.xlabel('Value')
plt.ylabel('Density')
plt.legend()
#plt.savefig('Probability density function'+'.pdf')
plt.savefig('Probability density function_Riverside.jpg',format='jpg', dpi=1000)
plt.show()


a0, loc0, scale0 = gamma.fit(BR, floc=0)  # We fix the location to zero
# print(a, loc, scale)
x0 = np.linspace(0, 1.5*max(BR), 1000)
pdf_fitted0 = gamma.pdf(x0, a0, loc0, scale0)
plt.figure(figsize=(10, 6))
plt.hist(BR, bins=3, density=True, color='blue', edgecolor='black', alpha=0.6)
plt.plot(x0, pdf_fitted0, 'r-', label=f'Gamma PDF\nshape={a:.2f}, scale={scale:.2f}')
plt.title('Histogram of Effective Biting Rates with Fitted Gamma PDF, Riverside County, California')
plt.xlabel('Normilized values of number of effetive Biting Rates')
plt.ylabel('Density')
plt.legend()
plt.savefig('Histogram_Bitingrate_Riverside_county_normalized.jpg',format='jpg', dpi=1000)
plt.show()








#------------------------------------------------------------------------------
#---------------------------plotting the estimated risk------------------------
#------------------------------------------------------------------------------
years = [year for year in range(2006, 2024)  if year != 2010]
for year in years:
    dama = list(get_temperature_for_year(file_path, year)[0]) 
    if year== 2021 or year==2022 or year==2023:
        if year==2021:
            Casetime=list(np.loadtxt("./RS2021.txt"))
            C=list(np.loadtxt("./C_esti_2021.txt"))
        elif year == 2022:
            Casetime=list(np.loadtxt("./RS2022.txt"))
            C=list(np.loadtxt("./C_esti_2022.txt"))
        elif year==2023:
            Casetime=list(np.loadtxt("./RS2023.txt"))
            C=list(np.loadtxt("./C_esti_2023.txt"))        
    else:
        Casetime=list(np.loadtxt(f"./RS{year}.txt"))
        dama = list(get_temperature_for_year(file_path, year)[0])
        C=list(np.loadtxt(f"./C{year}.txt"))   
    def nonzero_indices(lst):
        return [index for index, value in enumerate(lst) if value != 0]
    ct=nonzero_indices(Casetime)
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
        CB=1217596
        TH=x[0]+x[1]+x[2]+x[3]
        D=x[11]+x[12]+x[13]+x[14]
        ph=Coematrix[3,time]/(D+TH)
        dsHdt=(muHuman)-ph*Coematrix[7,time]*x[8]*x[0]-muHuman*x[0]
        #dsHdt=35-ph*Coematrix[7,time]*x[8]*x[0]-muHuman*x[0]
        deHdt=ph*Coematrix[7,time]*x[8]*x[0]-(incHumanWND)*x[1]-muHuman*x[1]
        diHdt=(incHumanWND)*x[1]-(infectionhumanWND)*x[2]-(muhumanWND)*x[2]-muHuman*x[2]
        drHdt=(infectionhumanWND)*x[2]-muHuman*x[3]
        #--------------------Mosquitoes' developement equations
        deggMdt=Coematrix[4,time]*(x[6]+x[7]+x[8])-Coematrix[5,time]*x[4]
        dacuMdt=Coematrix[5,time]*x[4]*max(0,(1-(x[5]/(C[time]))))-Coematrix[0,time]*Coematrix[1,time]*x[5]
        #-------------------Mosquitoe +pathogen+Bird+Human
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
#-----------------------------------------------------spill over test
    def find_first_nonzero(lst):
        # Iterate over the list and return the index of the first nonzero element
        for index, value in enumerate(lst):
            if value != 0:
                return index
        # Return -1 if no nonzero element is found
        return -1
    s=7*find_first_nonzero(Casetime)
    #s=[7*c for c in ct]
    #print("ct",ct)    
    x0=(2474000,0,0,0,1,1,1,1,0,100,100,100,100,100,100)
    r=2
    t=np.linspace(0,r*365,r*364)
    x=odeint(BMPPT,x0,t)
    R=x[(r-1)*365-1:r*365,3]
    I=x[(r-1)*365-1:r*365,2]
    M=x[(r-1)*365-1:r*365,6]
    m=M/max(M)


#----------------------------------------savinf result as txt (mosquitoe profiles)
    # df = pd.DataFrame(M, columns=['Numbers'])
    # df.to_csv(f"./Mosquitoe_Profile_Orange_County_{year}.txt", index=False, header=False, sep='\t')
#------------------------------------------------------------------------------


     
     
    for i in range(0,364):
        r=0 
        #plt.plot(i,M[i],Coematrix[3,i],'.k')
        #print(i,M[i],Coematrix[3,i],'.k')        
        A = contour_lines.collections[0].get_paths()
        for path in A:  # Loop through each path
            if path.contains_point((M[i], Coematrix[3,i])):  # Check if the point is within the path
                r=r+1
                break  # If the point is found in any path, you may want to break the loop
        A = contour_lines.collections[1].get_paths()
        for path in A:  # Loop through each path
            if path.contains_point((M[i], Coematrix[3,i])):  # Check if the point is within the path
                r=r+1
                break  # If the point is found in any path, you may want to break the loop
        A = contour_lines.collections[2].get_paths()
        for path in A:  # Loop through each path
            if path.contains_point((M[i], Coematrix[3,i])):  # Check if the point is within the path
                r=r+1
                break  # If the point is found in any path, you may want to break the loop
        if r==1:
            plt.axvline(i, color='y')
            # print('y')
        elif r==2:
            plt.axvline(i, color='orange') 
            # print('b')
        elif r==3:
            plt.axvline(i, color='r') 
            # print('r')
        elif r==0:
            plt.axvline(i, color='g')
            # print('g')              

    if year<2010:
        plt.title('Probability of spillover for 200'+str(f'{year}'))
        #plt.savefig('200' + str(f'{r}') + '.pdf')
        plt.savefig('200' + str(f'{year}'))
        plt.yticks(np.linspace(0, 1, 11))
        plt.plot(m, '--k', label='Estimated Profile of Mosquitoes')
        # for i in ct:
        #     if i == min(ct):
        #         plt.axvline(x=i*7, color='k', linestyle='--', label='Spillover Time')
        #     else:
        #         plt.axvline(x=i*7, color='k', linestyle='--')
        plt.axvline(s, color='b', linestyle='--', label='Spillover Time')
        plt.xlabel('Days')
        plt.ylabel('Estimated Number of Mosquitoes')
        plt.title('Estimated Risk - Riverside County - ' + str(year))
        plt.legend()
        plt.savefig(f'risk{year}.jpg', format='jpg', dpi=300)
        plt.show()
        #-------------------------------------------------------------------------
        file_name = f'm_profile_{year}.txt'
        with open(file_name, 'w') as file:
            for item in m:
                file.write(f"{item}\n")
        print(f"List saved to {file_name}")
        
    elif year == 2021 or year ==2022 or year==2023:
        plt.title('Probability of spillover for each day, 20'+str(f'{year}'))
        plt.yticks(np.linspace(0, 1, 11))
        plt.plot(m,'--k',label='estimated profile of mosquitoes')
        # for i in ct:
        #     if i==min(ct):
        #         plt.axvline(x=i*7, color='k', linestyle='--', label='Spillover time')
        #     else:
        #         plt.axvline(x=i*7, color='k', linestyle='--')
        plt.axvline(s, color='b', linestyle='--', label='spillover time')
        plt.xlabel('Days')
        plt.ylabel('Predicted number of mosquitoes')
        plt.title('Predicted Risk-Riverside County - '+str(year))
        plt.yticks(np.linspace(0, 1, 11))
        plt.legend()
        #plt.savefig(f'risk{year}.jpg', format='jpg', dpi=300)
        plt.show()
        #-------------------------------------------------------------------------
        file_name = f'm_profile_{year}.txt'
        with open(file_name, 'w') as file:
            for item in m:
                file.write(f"{item}\n")
        print(f"List saved to {file_name}")
        
    else:
        plt.title('Probability of spillover for each day, 20'+str(f'{year}'))
        plt.plot(m,'--k',label='estimated profile of mosquitoes')
        # for i in ct:
        #     if i==min(ct):
        #         plt.axvline(x=i*7, color='k', linestyle='--', label='Spillover time')
        #     else:
        #         plt.axvline(x=i*7, color='k', linestyle='--')
        plt.yticks(np.linspace(0, 1, 11))
        plt.axvline(s, color='b', linestyle='--', label='spillover time')
        plt.xlabel('Days')
        plt.ylabel('Estimated number of mosquitoes')
        plt.title('Estimated Risk-Riverside County - '+str(year))
        plt.yticks(np.linspace(0, 1, 11))
        plt.savefig(f'risk{year}.jpg', format='jpg', dpi=300)
        #plt.savefig(f'risk{year}'+'.pdf')
        plt.legend()
        plt.show()
        #-------------------------------------------------------------------------
        file_name = f'm_profile_{year}.txt'
        with open(file_name, 'w') as file:
            for item in m:
                file.write(f"{item}\n")
        print(f"List saved to {file_name}")
print("redto year =", redyear)
# print(redtogreen)
print("orange = ", oranyear)
print("yellow = ", yellyear)
# print(greenyear)
print("risk over year = ", risktoyear)
print("red to all = ",relativeredtorisk)
# plt.plot(redyear,'r--')
# plt.plot(oranyear,'o--')
# plt.plot(yellyear,'-y')
