# -*- coding: utf-8 -*-
"""
Created on Tue Jun 11 14:56:14 2024

@author: shosseini
"""

import numpy as np
from scipy.integrate import odeint
from matplotlib import pyplot as plt
from scipy.stats import gamma
import pandas as pd
from scipy.optimize import minimize
from scipy.stats import gamma, norm, multivariate_normal, chi2, kstest
import numpy as np
import matplotlib.pyplot as plt
from scipy.stats import norm, multivariate_normal, kstest
from copulas.visualization import scatter_2d
from copulas.multivariate import GaussianMultivariate
from scipy.stats import linregress, kstest


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
#C2000 = list(np.loadtxt("./C2000.txt"))
T2000 = list(np.loadtxt("./T2000.txt"))
H2000 = list(np.loadtxt("./H2000.txt"))

#C2001 = list(np.loadtxt("./C2001.txt"))
T2001 = list(np.loadtxt("./T2001.txt"))
H2001 = list(np.loadtxt("./H2001.txt"))

#C2002 = list(np.loadtxt("./C2002.txt"))
T2002 = list(np.loadtxt("./T2002.txt"))
H2002 = list(np.loadtxt("./H2002.txt"))

#C2003 = list(np.loadtxt("./C2003.txt"))
T2003 = list(np.loadtxt("./T2003.txt"))
H2003 = list(np.loadtxt("./H2003.txt"))

#C2004 = list(np.loadtxt("./C2004.txt"))
T2004 = list(np.loadtxt("./T2004.txt"))
H2004 = list(np.loadtxt("./H2004.txt"))

#C2005 = list(np.loadtxt("./C2005.txt"))
T2005 = list(np.loadtxt("./T2005.txt"))
H2005 = list(np.loadtxt("./H2005.txt"))






#--------------------------------------Section of importing data 
C2006 = list(np.loadtxt("./Cf2006.txt"))
T2006 = list(np.loadtxt("./T2006.txt"))
H2006 = list(np.loadtxt("./H2006.txt"))
#---------------------------------------------------------------
C2007 = list(np.loadtxt("./Cf2007.txt"))
T2007 = list(np.loadtxt("./T2007.txt"))
H2007 = list(np.loadtxt("./H2007.txt"))
#---------------------------------------------------------------
C2008 = list(np.loadtxt("./Cf2008.txt"))
T2008 = list(np.loadtxt("./T2008.txt"))
H2008 = list(np.loadtxt("./H2008.txt"))
#---------------------------------------------------------------
C2009 = list(np.loadtxt("./Cf2009.txt"))
T2009 = list(np.loadtxt("./T2009.txt"))
H2009 = list(np.loadtxt("./H2009.txt"))
#---------------------------------------------------------------
C2011 = list(np.loadtxt("./Cf2011.txt"))
T2011 = list(np.loadtxt("./T2011.txt"))
H2011 = list(np.loadtxt("./H2011.txt"))
#---------------------------------------------------------------
C2012 = list(np.loadtxt("./Cf2012.txt"))
T2012 = list(np.loadtxt("./T2012.txt"))
H2012 = list(np.loadtxt("./H2012.txt"))
#---------------------------------------------------------------
C2013 = list(np.loadtxt("./Cf2013.txt"))
T2013 = list(np.loadtxt("./T2013.txt"))
H2013 = list(np.loadtxt("./H2013.txt"))
#---------------------------------------------------------------
C2014 = list(np.loadtxt("./Cf2014.txt"))
T2014 = list(np.loadtxt("./T2014.txt"))
H2014 = list(np.loadtxt("./H2014.txt"))
#---------------------------------------------------------------
C2015 = list(np.loadtxt("./Cf2015.txt"))
T2015 = list(np.loadtxt("./T2015.txt"))
H2015 = list(np.loadtxt("./H2015.txt"))
#---------------------------------------------------------------
C2016 = list(np.loadtxt("./Cf2016.txt"))
T2016 = list(np.loadtxt("./T2016.txt"))
H2016 = list(np.loadtxt("./H2016.txt"))
#---------------------------------------------------------------
C2017 = list(np.loadtxt("./Cf2017.txt"))
T2017 = list(np.loadtxt("./T2017.txt"))
H2017 = list(np.loadtxt("./H2017.txt"))
#---------------------------------------------------------------
C2018 = list(np.loadtxt("./Cf2018.txt"))
T2018 = list(np.loadtxt("./T2018.txt"))
H2018 = list(np.loadtxt("./H2018.txt"))
#---------------------------------------------------------------
C2019 = list(np.loadtxt("./Cf2019.txt"))
T2019 = list(np.loadtxt("./T2019.txt"))
H2019 = list(np.loadtxt("./H2019.txt"))
#---------------------------------------------------------------
C2020 = list(np.loadtxt("./Cf2020.txt"))
T2020 = list(np.loadtxt("./T2020.txt"))
H2020 = list(np.loadtxt("./H2020.txt"))
#---------------------------------------------------------------
T2021 = list(np.loadtxt("./T2021.txt"))
H2021 = list(np.loadtxt("./H2021.txt"))
Cf2021 = list(np.loadtxt("./C_pre_2021.txt"))
#---------------------------------------------------------------
#---------------------------------------------------------------
T2022 = list(np.loadtxt("./T2022.txt"))
H2022 = list(np.loadtxt("./H2022.txt"))
Cf2022 = list(np.loadtxt("./C_pre_2022.txt"))
#---------------------------------------------------------------
T2023 = list(np.loadtxt("./T2023.txt"))
H2023 = list(np.loadtxt("./H2023.txt"))
Cf2023 = list(np.loadtxt("./C_pre_2023.txt"))
#--------------------------------------Section of the coeficient of
#                               the nonlinear model for pipiens mosquitoes----------------


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


# we will fill this matrix by coeficient of each day so we have 365 set of data
#------------------------------------------------------------------------------
Coematrix=np.zeros([9,365])    
N=[]
BR=[]
B=[]



# #-----------------------------ploting profile of the mosquitoe for each year---
# #-------------Also it collects the sample of number of effective populations---
# #------------------------------------------------------------------------------

# years = [year for year in range(2006, 2021) if year != 2010]
# for year in years:
#     redyear=[]
#     oranyear=[]
#     yellyear=[]
#     redtogreen=[]
#     greenyear=[]
#     risktoyear=[]
#     relativeredtorisk=[]
#     dama = list(np.loadtxt(f"./T{year}.txt"))
#     C=list(np.loadtxt(f"./Cf{year}.txt"))
#     Casetime=list(np.loadtxt(f"./orangecases{year}.txt"))
#     def nonzero_indices(lst):
#         return [index for index, value in enumerate(lst) if value != 0]
#     ct=nonzero_indices(Casetime)
#     #print(ct)
#     for time in range(0,365):
#         Coematrix[0,time]=gammaMP(time)
#         Coematrix[1,time]=pLSP(time)
#         Coematrix[2,time]=muMP(time)
#         Coematrix[3,time]=BRP(time)
#         Coematrix[4,time]=alphaEP(time)
#         Coematrix[5,time]=bEP(time)
#         Coematrix[6,time]=PDRP(time)
#         Coematrix[7,time]=VCP(time)

#     def BMPPT(x,time):
#         time=int(divmod(time, 365)[1])
#         x[0]=max(0,x[0]);x[1]=max(0,x[1]);x[2]=max(0,x[2]);x[3]=max(0,x[3]);x[4]=max(0,x[4])
#         x[5]=max(0,x[5]);x[6]=max(0,x[6]);x[7]=max(0,x[7]);x[8]=max(0,x[8]);x[9]=max(0,x[9]);
#         x[11]=max(0,x[11]);x[12]=max(0,x[12]);x[13]=max(0,x[13]);x[14]=max(0,x[14])
#         #**************************************************************************
#         sigma=15;mu=110
#         #FB=1.96*(1/np.sqrt((2*np.pi)*(sigma**2)))*np.exp(-(2*sigma**(-2))*(((np.divmod(t,365)[1])-mu)**2))
#         FB=2
#         EVB=0.021
#         muE=4.9315*(10**(-4))
#         BDR=9.041*(10**(-4))
#         muF=6.30136*(10**(-4))
#         muB=(0.15/365)
#         incB=1/3
#         recoveryrateC=1/6
#         muwnvB=0.9
#         muHbirth=5
#         muhumanWND=0.01
#         incHumanWND=1/6
#         infectionhumanWND=1/3
#         muHuman=(1/(77*365))
#         #***** equations related to dynamic of Culex pipiens***********************
#         #kC=D
#         #--------------------Human equations
#         CB=20000
#         TH=x[0]+x[1]+x[2]+x[3]
#         D=x[11]+x[12]+x[13]+x[14]
#         ph=Coematrix[3,time]/(D+TH)
#         dsHdt=(muHuman*(TH))-ph*Coematrix[7,time]*x[8]*x[0]-muHuman*x[0]
#         #dsHdt=35-ph*Coematrix[7,time]*x[8]*x[0]-muHuman*x[0]
#         deHdt=ph*Coematrix[7,time]*x[8]*x[0]-(incHumanWND)*x[1]-muHuman*x[1]
#         diHdt=(incHumanWND)*x[1]-(infectionhumanWND)*x[2]-(muhumanWND)*x[2]-muHuman*x[2]
#         drHdt=(infectionhumanWND)*x[2]-muHuman*x[3]
#         #--------------------Mosquitoes' developement equations
#         deggMdt=Coematrix[4,time]*(x[6]+x[7]+x[8])-Coematrix[5,time]*x[4]
#         dacuMdt=Coematrix[5,time]*x[4]*max(0,(1-(x[5]/(C[time]))))-Coematrix[0,time]*Coematrix[1,time]*x[5]
#         #-------------------Mosquitoe +pathogen+Bird+Human
#         dsMdt=Coematrix[0,time]*Coematrix[1,time]*x[5]-ph*x[13]*x[6]-Coematrix[2,time]*x[6]
#         deMdt=ph*x[13]*x[6]-Coematrix[6,time]*x[7]-Coematrix[2,time]*x[7]
#         diMdt=Coematrix[6,time]*x[7]-Coematrix[2,time]*x[8]
#         #------------------Birds' developement equations
#         deggBdt=FB*(x[11]+x[12]+x[13]+x[14])-EVB*x[9]-muE*x[9]
#         dfleBdt=EVB*x[9]*max(0,(1-(x[10]/CB)))-BDR*x[10]-muF*x[10]
#         #------------------Bird+pathogen+Mosquito
#         dsBdt=BDR*x[10]-ph*Coematrix[7,time]*x[8]*x[11]-muB*x[11]
#         deBdt=ph*Coematrix[7,time]*x[8]*x[11]-incB*x[12]-muB*x[12]
#         diBdt=incB*x[12]-muwnvB*x[13]-muB*x[13]-recoveryrateC*x[13]
#         drBdt=recoveryrateC*x[13]-muB*x[14]
#         dxdt=(dsHdt,deHdt,diHdt,drHdt,deggMdt,dacuMdt,dsMdt,deMdt,diMdt,deggBdt,dfleBdt,dsBdt,deBdt,diBdt,drBdt)
#         return(dxdt)
#     t=np.linspace(0,365,364)
#     x0=(3000000,0,0,0,10,10,10,10,10,100,100,100,100,100,100)
#     x=odeint(BMPPT,x0,t)
#     muB=(0.15/365)
#     incB=1/3
#     rC=1/6
#     muwnvB=0.9
#     R=x[0:365,3]
#     I=x[0:365,2]
#     M=x[0:365,6]
#     BS=x[0:365,11]
#     TH=x[0:365,0]+x[0:365,1]+x[0:365,2]+x[0:365,3]
#     D=x[0:365,11]+x[0:365,12]+x[0:365,13]+x[0:365,14]
#     BetaBM=Coematrix[7,0:364]*Coematrix[3,0:364]/(D+TH)
#     BetaMB=Coematrix[3,0:364]/(D+TH)
#     R2=(BetaBM*BetaMB*M*BS*Coematrix[6,0:364]*incB)/((muB+incB)*(rC+muwnvB+muB)*(Coematrix[2,0:364])*(Coematrix[6,0:364]+Coematrix[2,0:364]))
#     R0=[np.sqrt(R2[i]) for i in range(0,364)]

   
#     fig, axes = plt.subplots(1, 2, figsize=(14, 6))

#     # Plot 1: Estimated profile of mosquitoes
#     axes[0].plot(M / max(M), '--k', label='Estimated profile of mosquitoes')
#     for i in ct:
#         if i == min(ct):
#             axes[0].axvline(x=i * 7, color='r', linestyle='-', label='Spillover time')
#         else:
#             axes[0].axvline(x=i * 7, color='r', linestyle='--')
#     axes[0].set_xlabel('Days')
#     axes[0].set_ylabel('Estimated relative abundance of mosquitoes')
#     axes[0].set_title('Orange County, California - ' + str(year))
#     axes[0].legend()
#     axes[0].grid(True)
    
#     # Plot 2: Basic reproductive number
#     axes[1].plot(R0, '--k', label='Basic reproductive number')
#     for i in ct:
#         if i == min(ct):
#             axes[1].axvline(x=i * 7, color='r', linestyle='-', label='Spillover time')
#         else:
#             axes[1].axvline(x=i * 7, color='r', linestyle='--')
#     axes[1].set_xlabel('Days')
#     axes[1].set_ylabel('Estimated basic reproductive number')
#     axes[1].set_title('Orange County, California - ' + str(year))
#     axes[1].legend()
#     axes[1].grid(True)
    
#     plt.tight_layout()
#     plt.savefig(f'Profile_and_R0_orange_county_{year}.pdf', format='pdf', dpi=1000)
#     plt.show()
    
#     BR.append(Coematrix[3,ct[0]*7])
#     N.append(M[ct[0]*7])
#     B.append(R0[ct[0]*7])
#     # for i in ct:
#     # #     print(M[i*7])
#     #     B.append(R0[i*7])
#     #     N.append(M[i*7])

# #---------------------------------------------Copulat joint distribution-------

# print("BR_array =", np.array(B))
# print("N_array =",np.array(N))
# # Provided data
# BR_array = np.array(B)
# N_array =np.array(N)
'''
BR_array = [1.53293217 ,2.82375948, 3.03207687, 2.45577936, 1.71378965, 2.74416761,
 2.10249095, 2.69200256, 2.55298161, 2.62414203, 1.92365222, 1.64652918,
 3.66656926, 3.03662025]
N_array = [1.69628696e+10, 2.69182731e+10, 4.51135413e+10, 2.56554904e+10,
 4.92227175e+10, 4.14478518e+10, 2.49309124e+10, 4.45614678e+10,
 2.57437863e+10, 2.55246120e+10 ,2.81278301e+10 ,3.20573836e+10,
 4.05848078e+10 ,3.24250238e+10]

# Fit Gaussian distribution to data
def fit_gaussian(data):
    mean, std = norm.fit(data)
    return mean, std

mean_N, std_N = fit_gaussian(N_array)
mean_BR, std_BR = fit_gaussian(BR_array)

# Transform data using the fitted normal CDF
F_N = norm.cdf(N_array, mean_N, std_N)
F_BR = norm.cdf(BR_array, mean_BR, std_BR)

# Apply probability integral transform (PIT)
U = norm.ppf(F_N)
V = norm.ppf(F_BR)

# Calculate correlation coefficient and fit Gaussian copula
correlation_coefficient = np.corrcoef(BR_array, N_array)[0, 1]
mean = [0, 0]
cov = [[1, correlation_coefficient], [correlation_coefficient, 1]]
rv = multivariate_normal(mean, cov)

# Generate grid for plotting
n = np.linspace(0, 10 * max(N_array), 1000)
br = np.linspace(0, 10 * max(BR_array), 1000)
N_grid, BR_grid = np.meshgrid(n, br)

U_grid = norm.ppf(norm.cdf(N_grid, mean_N, std_N))
V_grid = norm.ppf(norm.cdf(BR_grid, mean_BR, std_BR))
pos = np.dstack((U_grid, V_grid))
Z = rv.pdf(pos)
probability_levels = [0.8, 0.9, 0.95]
chi2_values = [chi2.ppf(level, df=2) for level in probability_levels]
contour_levels = [rv.pdf([0, 0]) * np.exp(-0.5 * chi2_val) for chi2_val in chi2_values]
fig, ax = plt.subplots(figsize=(12, 8))
contour_filled = ax.contourf(N_grid, BR_grid, Z, levels=40, cmap='viridis')
plt.colorbar(contour_filled, ax=ax, label='Joint PDF')
contour_lines = plt.contour(N_grid, BR_grid, Z, levels=sorted(contour_levels), colors=['black', 'black', 'black'])
ax.clabel(contour_lines,  fmt={contour_levels[0]: probability_levels[0], contour_levels[1]: probability_levels[1], contour_levels[2]: probability_levels[2]}, inline=True, fontsize=10)


ax.set_xlabel('Effective Value of Mosquito Population')
ax.set_ylabel('Effective Value of Basic Reproductive Number')
ax.set_title('Joint PDF')
plt.grid(True)
plt.xlim(0, 1.3 * max(N_array))
plt.ylim(0, 1.3 * max(BR_array))
plt.savefig('Joint Copula pdf Orange County' + '.pdf')
plt.show()

# Plot histograms with fitted normal distributions
fig, axes = plt.subplots(1, 2, figsize=(14, 6))

# Histogram for BR_array with Gaussian fit
axes[0].hist(BR_array, bins=6, density=True, alpha=0.6, color='g', label='BR_array')
xmin, xmax = axes[0].get_xlim()
x = np.linspace(xmin, xmax, 100)
p = norm.pdf(x, mean_BR, std_BR)
axes[0].plot(x, p, 'k', linewidth=2)
title = f'Gaussian Fit: mu = {mean_BR:.2f}, std = {std_BR:.2f}'
axes[0].set_title(title)
axes[0].set_xlabel('BR_array')
axes[0].set_ylabel('Density')
axes[0].legend()

# Histogram for N_array with Gaussian fit
axes[1].hist(N_array, bins=6, density=True, alpha=0.6, color='b', label='N_array')
xmin, xmax = axes[1].get_xlim()
x = np.linspace(xmin, xmax, 100)
p = norm.pdf(x, mean_N, std_N)
axes[1].plot(x, p, 'k', linewidth=2)
title = f'Gaussian Fit: mu = {mean_N:.2e}, std = {std_N:.2e}'
axes[1].set_title(title)
axes[1].set_xlabel('N_array')
axes[1].set_ylabel('Density')
axes[1].legend()

plt.tight_layout()
plt.savefig('Histograms with Fitted Normal Distributions' + '.pdf')
plt.show()

# Goodness-of-fit tests
# Kolmogorov-Smirnov test for marginals
ks_stat_N, p_value_N = kstest(N_array, 'norm', args=(mean_N, std_N))
ks_stat_BR, p_value_BR = kstest(BR_array, 'norm', args=(mean_BR, std_BR))

print(f'KS test for N_array: stat={ks_stat_N}, p-value={p_value_N}')
print(f'KS test for BR_array: stat={ks_stat_BR}, p-value={p_value_BR}')

# Cramér-von Mises statistic for empirical vs. fitted copula
empirical_copula_values = rv.cdf(np.column_stack((F_N, F_BR)))
theoretical_copula_values = rv.cdf(np.column_stack((U, V)))

# Ensure both are 2-dimensional arrays for element-wise operations
cvm_stat = np.sum((empirical_copula_values - theoretical_copula_values) ** 2)
print(f'Cramer-von Mises statistic: {cvm_stat}')






#------------------------------------------------------------------------------
#---------------------------plotting the estimated risk------------------------
#------------------------------------------------------------------------------
file_path = 'info.csv'
RRR=[]
EEE=[]
YYY=[]
years = [year for year in range(1981, 2024)]
#years = [year for year in range(2006, 2024)]
for year in years:
    if year>=1981 and year<=2000:
        Casetime=list(np.loadtxt(f"./orangecases{2006}.txt"))
        C=list(np.loadtxt(f"./C_pre_{year}.txt"))
        d0=get_temperature_for_year(file_path, year)[0]
        dama=[d0[i] for i in range(0,365)]
    elif year>=2000 and year<=2005:
        Casetime=list(np.loadtxt(f"./orangecases{2006}.txt"))
        C=list(np.loadtxt(f"./C_pre_{year}.txt"))
        dama = list(np.loadtxt(f"./T{year}.txt"))
    elif year==2021:
        Casetime=list(np.loadtxt(f"./orangecases{year}.txt"))
        C=list(np.loadtxt("./C_pre_2021.txt"))
        dama = list(np.loadtxt(f"./T{year}.txt"))
    elif year==2022:
        Casetime=list(np.loadtxt(f"./orangecases{year}.txt"))
        C=list(np.loadtxt("./C_pre_2022.txt"))
        dama = list(np.loadtxt(f"./T{year}.txt"))
    elif year==2023:
        Casetime=list(np.loadtxt(f"./orangecases{year}.txt"))
        C=list(np.loadtxt("./C_pre_2023.txt"))
        dama = list(np.loadtxt(f"./T{year}.txt"))
    else:
        Casetime=list(np.loadtxt(f"./orangecases{year}.txt"))
        dama = list(np.loadtxt(f"./T{year}.txt"))
        #C=list(np.loadtxt(f"./Cf{year}.txt"))
        C=list(np.loadtxt(f"./C_pre_{year}.txt"))
    
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
#-----------------------------------------------------spill over test
    def find_first_nonzero(lst):
        # Iterate over the list and return the index of the first nonzero element
        for index, value in enumerate(lst):
            if value != 0:
                return index
        # Return -1 if no nonzero element is found
        return -1
    s=7*find_first_nonzero(Casetime)
    #s=7*Casetime
    #print(s)
    # t=np.linspace(0,365,364)
    # x0=(3000000,0,0,0,10,10,10,10,10,100,100,100,100,100,100)
    # x=odeint(BMPPT,x0,t)
    # R=x[0:365,3]
    # I=x[0:365,2]
    # M=x[0:365,6]
    
    t=np.linspace(0,365,364)
    x0=(3000000,0,0,0,10,10,10,10,10,100,100,100,100,100,100)
    x=odeint(BMPPT,x0,t)
    muB=(0.15/365)
    incB=1/3
    rC=1/6
    muwnvB=0.9
    R=x[0:365,3]
    I=x[0:365,2]
    M=x[0:365,6]
    BS=x[0:365,11]
    TH=x[0:365,0]+x[0:365,1]+x[0:365,2]+x[0:365,3]
    D=x[0:365,11]+x[0:365,12]+x[0:365,13]+x[0:365,14]
    BetaBM=Coematrix[7,0:364]*Coematrix[3,0:364]/(D+TH)
    BetaMB=Coematrix[3,0:364]/(D+TH)
    R2=(BetaBM*BetaMB*M*BS*Coematrix[6,0:364]*incB)/((muB+incB)*(rC+muwnvB+muB)*(Coematrix[2,0:364])*(Coematrix[6,0:364]+Coematrix[2,0:364]))
    R0=[np.sqrt(R2[i]) for i in range(0,364)]
#----------------------------------------savinf result as txt (mosquitoe profiles)
    # df = pd.DataFrame(M, columns=['Numbers'])
    # df.to_csv(f"./Mosquitoe_Profile_Orange_County_{year}.txt", index=False, header=False, sep='\t')
#------------------------------------------------------------------------------
#--------------------------------------------------------------------------------------
    yellow=0
    green=0
    orange=0
    red=0
    
    for i in range(0,364):
        r=0 
        #plt.plot(i,M[i],Coematrix[3,i],'.k')
        #print(i,M[i],Coematrix[3,i],'.k')        
        A = contour_lines.collections[0].get_paths()
        for path in A:  # Loop through each path
            if path.contains_point((M[i], R0[i])):  # Check if the point is within the path
                r=r+1
                break  # If the point is found in any path, you may want to break the loop
        A = contour_lines.collections[1].get_paths()
        for path in A:  # Loop through each path
            if path.contains_point((M[i], R0[i])):  # Check if the point is within the path
                r=r+1
                break  # If the point is found in any path, you may want to break the loop
        A = contour_lines.collections[2].get_paths()
        for path in A:  # Loop through each path
            if path.contains_point((M[i], R0[i])):  # Check if the point is within the path
                r=r+1
                break  # If the point is found in any path, you may want to break the loop
        if r==1:
            plt.axvline(i, color='y')
            yellow=yellow+1
        elif r==2:
            plt.axvline(i, color='orange') 
            orange=orange+1
        elif r==3:
            plt.axvline(i, color='r') 
            red=red+1
        elif r==0:
            plt.axvline(i, color='g')
            green=green+1
    RRR.append(red)
    EEE.append(red/(yellow+orange+red))
    YYY.append(yellow+orange+red)
            
    if year>2020:
        plt.plot(M/max(M),'--k',label='Predicted profile of mosquitoes')
        plt.axvline(s, color='b', linestyle='--', label='Spillover time')
        plt.xlabel('Days')
        plt.ylabel('Predicted realtive abundance of mosquitoes')
        plt.title('Predicted risk-Orange County, California - '+str(year))
        #plt.savefig(f'risk{year}.jpg', format='jpg', dpi=300)
        plt.legend()
        plt.savefig(f'predicted_risk{year}'+'.pdf')
        plt.show()
        
    else:
         #plt.legend(['green: P < 0.03', 'yellow: 0.03 < P < 0.3', 'orange: 0.3 < P < 0.6', 'red: 0.6 < P '],loc='upper left')
         # plt.yticks([])
        plt.plot(M/max(M),'--k',label='Estimated profile of mosquitoes')
        plt.axvline(s, color='b', linestyle='--', label='Spillover time')
        plt.xlabel('Days')
        plt.ylabel('Estimated realtive abundance of mosquitoes')
        plt.title('Estimated risk-Orange County, California - '+str(year))
        #plt.savefig(f'risk{year}.jpg', format='jpg', dpi=300)
        plt.legend()
        plt.savefig(f'risk{year}'+'.pdf')
        plt.show()
    #--------------------------------------------------------------------
    # file_name = f'm_profile_{year}.txt'
    # with open(file_name, 'w') as file:
    #     for item in M:
    #         file.write(f"{item}\n")
    # print(f"List saved to {file_name}")
'''
redtoyear = [110, 81, 98, 107, 91, 84, 79, 103, 92, 109, 69, 87, 69, 90, 96, 95, 87, 80, 84, 77, 97, 79, 99, 76, 77, 105, 99, 101, 88, 83, 85, 97, 101, 109, 104, 96, 109, 96, 93, 98, 109, 114, 90]
riskoveryear =  [132, 100, 120, 138, 118, 113, 125, 125, 119, 141, 113, 123, 115, 114, 120, 131, 137, 95, 123, 122, 129, 114, 120, 118, 107, 130, 119, 136, 115, 108, 107, 123, 124, 156, 148, 132, 136, 132, 117, 140, 140, 139, 110]
redtorisk =  [0.8333333333333334, 0.81, 0.8166666666666667, 0.7753623188405797, 0.7711864406779662, 0.7433628318584071, 0.632, 0.824, 0.773109243697479, 0.7730496453900709, 0.6106194690265486, 0.7073170731707317, 0.6, 0.7894736842105263, 0.8, 0.7251908396946565, 0.635036496350365, 0.8421052631578947, 0.6829268292682927, 0.6311475409836066, 0.751937984496124, 0.6929824561403509, 0.825, 0.6440677966101694, 0.719626168224299, 0.8076923076923077, 0.8319327731092437, 0.7426470588235294, 0.7652173913043478, 0.7685185185185185, 0.794392523364486, 0.7886178861788617, 0.8145161290322581, 0.6987179487179487, 0.7027027027027027, 0.7272727272727273, 0.8014705882352942, 0.7272727272727273, 0.7948717948717948, 0.7, 0.7785714285714286, 0.8201438848920863, 0.8181818181818182]
redtoyear = [a/365 for a in (redtoyear[10:len(redtoyear)])]
redtorisk = [a/365 for a in (redtorisk[10:len(redtorisk)])]
riskoveryear=[a/365 for a in (riskoveryear[10:len(riskoveryear)])]
#-------------------------------------------------------------------------------
x = np.arange(len(range(1991, 2024)))  # Indices corresponding to the years from 1989 to 2023
y = redtoyear

# Perform linear regression
slope, intercept, r_value, p_value, std_err = linregress(x, y)
regression_line = slope * x + intercept
years = np.arange(1991, 2024)  # Years from 1989 to 2023
slope=np.round(slope,decimals=5)
plt.figure(figsize=(12, 6))
plt.plot(x, y, '--o', label='Annual high-risk indicator')
plt.plot(x, regression_line, 'r', label=f'Regression line: y = {slope:f}x + {intercept:f}')
# plt.xlabel('Year')
plt.ylabel('Annual high-risk indicator')
plt.title('Linear trend of annual high-risk indicator over the years, Orange County, California')
plt.legend()
plt.xticks(ticks=x, labels=years, rotation=45)
plt.savefig('The annual high-risk indicator over the years'+'.pdf')
plt.show()
#-------------------------------------------------------------------------------------
x = np.arange(len(range(1991, 2024)))  # Indices corresponding to the years from 1989 to 2023
y = redtorisk
slope, intercept, r_value, p_value, std_err = linregress(x, y)
regression_line = slope * x + intercept
years = np.arange(1991, 2024)  # Years from 1989 to 2023
slope=np.round(slope,decimals=5)
plt.figure(figsize=(12, 6))
plt.plot(x, y, '--o', label='Annual relative high-risk indicator')
plt.plot(x, regression_line, 'r', label=f'Regression line: y = {slope:f}x + {intercept:f}')
# plt.xlabel('Year')
plt.ylabel('annual relative high-risk indicator')
plt.title('Linear trend of annual relative high-risk indicator over the years, Orange County, California')
plt.legend()
plt.xticks(ticks=x, labels=years, rotation=45)
plt.savefig('The annual relative high-risk indicator'+'.pdf')
plt.show()


def plot_regression_and_test(x, y, title):
    slope, intercept, r_value, p_value, std_err = linregress(x, y)
    regression_line = slope * np.array(x) + intercept
    residuals = np.array(y) - regression_line
    ks_stat, ks_p_value = kstest(residuals, 'norm', args=(np.mean(residuals), np.std(residuals)))
    print(f'{title} - Linear Regression:')
    print(f'Slope: {slope:f}')
    print(f'Intercept: {intercept:f}')
    print(f'R-squared: {r_value**2:f}')
    print(f'p-value: {p_value:f}')
    print(f'Std Err: {std_err:f}')
    print(f'KS test stat: {ks_stat:f}, KS test p-value: {ks_p_value:f}')
    print()

# Plot regression and test for each time series
plot_regression_and_test(list(range(len(redtoyear))), redtoyear, 'Red to Year')
#plot_regression_and_test(list(range(len(riskoveryear))), riskoveryear, 'Risk over Year')
plot_regression_and_test(list(range(len(redtorisk))), redtorisk, 'Red to Risk')

