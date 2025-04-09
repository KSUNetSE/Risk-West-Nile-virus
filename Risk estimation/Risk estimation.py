# -*- coding: utf-8 -*-
"""
Created on Tue Apr 16 11:26:55 2024

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

def moving_average(vector, window_size):
    return np.convolve(vector, np.ones(window_size)/window_size, mode='valid')


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
T2021 = list(np.loadtxt("./2021forecasted_tem.txt"))
H2021 = list(np.loadtxt("./2021forecasted_hum.txt"))
Cf2021 = list(np.loadtxt("./C_esti_2021.txt"))
#---------------------------------------------------------------
#---------------------------------------------------------------
T2022 = list(np.loadtxt("./T2022.txt"))
H2022 = list(np.loadtxt("./H2022.txt"))
Cf2022 = list(np.loadtxt("./C_esti_2022.txt"))
#---------------------------------------------------------------
T2023 = list(np.loadtxt("./T2023.txt"))
H2023 = list(np.loadtxt("./H2023.txt"))
Cf2023 = list(np.loadtxt("./C_esti_2023.txt"))
#--------------------------------------Section of the coeficient of
#                               the nonlinear model for pipiens mosquitoes----------------




years = [2006, 2007, 2008, 2009, 2011, 2012, 2013, 2014, 2015, 2016, 2017, 2018, 2019, 2020]

# Load all C{year} files into a 2D NumPy array
C_values = np.array([np.loadtxt(f"./Cf{year}.txt") for year in years])

# Compute the element-wise mean across all years
meanC = np.mean(C_values, axis=0)

# Convert to a list if needed
meanC = meanC.tolist()

print("mean",meanC)







def tim2tem(time):
    tempc=dama[time]
    #if tempc>25:tempc=25
    # tim2tem=(t2t-32)*(5/9)
    return(tempc)
# def tim2tem(time):
#     if time<180:
#         tempc=dama[time]-7
#     else:
#         tempc=dama[time]-4  
#     # if tempc>25:tempc=25
#     # tim2tem=(t2t-32)*(5/9)
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



#-----------------------------ploting profile of the mosquitoe for each year---
#-------------Also it collects the sample of number of effective populations---
#------------------------------------------------------------------------------

years = [year for year in range(2006, 2021) if year != 2010]
for year in years:
    redyear=[]
    oranyear=[]
    yellyear=[]
    redtogreen=[]
    greenyear=[]
    risktoyear=[]
    relativeredtorisk=[]
    dama = list(np.loadtxt(f"./T{year}.txt"))
    C=list(np.loadtxt(f"./Cf{year}.txt"))
    Casetime=list(np.loadtxt(f"./orangecases{year}.txt"))

    





    def get_first_nonzero_and_repeat(lst):
        """
        Finds the first non-zero value in a list, returns a list repeating its index.
    
        Parameters:
        lst (list): List of numbers.
    
        Returns:
        list: [index] repeated 'value' times.
        """
        for i, value in enumerate(lst):
            if value != 0:
                return [i] * int(value)  # Convert to int
        return []  # Return empty list if all values are zero
    def max_nonzero_index(vec):
        """
        Returns the index of the largest nonzero value in the vector.
        If all values are the same, returns the first index.
        
        Parameters:
            vec (list or np.ndarray): Input vector.
        
        Returns:
            int: Index of the largest nonzero value (or first index if all values are the same).
        """
        vec = np.array(vec)  # Ensure it's a NumPy array
        
        # Check if all values are the same
        if np.all(vec == vec[0]):  
            return 0  # Return the first index
        
        # Find indices of nonzero elements
        nonzero_indices = np.nonzero(vec)[0]
        
        if nonzero_indices.size == 0:  # If all elements are zero, return first index
            return 0  
        
        # Extract nonzero values
        nonzero_values = vec[nonzero_indices]
        
        # Find the index of the largest nonzero value
        max_index = nonzero_indices[np.argmax(nonzero_values)]
        
        return int(max_index )
    ct=get_first_nonzero_and_repeat(Casetime)
    ctmax=max_nonzero_index(Casetime)
    print("ct",ct)
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
    t=np.linspace(0,365,364)
    x0=(3000000,0,0,0,10,10,10,10,10,100,100,100,100,100,100)
    x=odeint(BMPPT,x0,t)
    muB=(0.15/365)
    incB=1/3
    rC=1/6
    muwnvB=0.9
    R=x[0:365,3]
    I=x[0:365,2]
    M=np.cumsum(x[0:365,6])
    Mos=x[0:365,6]
    BS=x[0:365,11]
    TH=x[0:365,0]+x[0:365,1]+x[0:365,2]+x[0:365,3]
    D=x[0:365,11]+x[0:365,12]+x[0:365,13]+x[0:365,14]
    BetaBM=Coematrix[7,0:364]*Coematrix[3,0:364]/(D+TH)
    BetaMB=Coematrix[3,0:364]/(D+TH)
    R2=(BetaBM*BetaMB*Mos*BS*Coematrix[6,0:364]*incB)/((muB+incB)*(rC+muwnvB+muB)*(Coematrix[2,0:364])*(Coematrix[6,0:364]+Coematrix[2,0:364]))
    R0=[np.sqrt(R2[i]) for i in range(0,364)]
    basicR0=R0
    
    
    
    R0=np.cumsum(R0)

   
    fig, axes = plt.subplots(1, 2, figsize=(14, 6))

    # Plot 1: Estimated profile of mosquitoes
    axes[0].plot(Mos / max(Mos), '--k', label='Estimated profile of mosquitoes')


    axes[0].axvline(x=ctmax*7, color='r', linestyle='-', label='Spillover time')

    axes[0].set_xlabel('Days')
    axes[0].set_ylabel('Estimated relative abundance of mosquitoes')
    axes[0].set_title('Orange County, California - ' + str(year))
    axes[0].legend()
    axes[0].grid(True)
    
    #Plot 2: Basic reproductive number
    axes[1].plot(basicR0, '--k', label='Basic reproductive number')
    
    axes[1].axvline(x=(ctmax) * 7, color='r', linestyle='-', label='Spillover time')
    
    axes[1].set_xlabel('Days')
    axes[1].set_ylabel('Estimated basic reproductive number')
    axes[1].set_title('Orange County, California - ' + str(year))
    axes[1].legend()
    axes[1].grid(True)
    
    plt.tight_layout()
    #plt.savefig(f'Profile_Orange_county{year}.jpg', format='jpg', dpi=1000)
    plt.show()
    
    
    
    
    # BR.append(Coematrix[3,ct[0]*7])
    # N.append(M[ct[0]*7])
    # B.append(R0[ct[0]*7])
    for i in ct:
        #print(M[i*7])
        B.append(R0[i*7])
        N.append(M[i*7])

#---------------------------------------------Copulat joint distribution-------


# Provided data
BR_array = np.array(B)
N_array =np.array(N)

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
probability_levels = [0.85, 0.9, 0.95]
chi2_values = [chi2.ppf(level, df=2) for level in probability_levels]
contour_levels = [rv.pdf([0, 0]) * np.exp(-0.5 * chi2_val) for chi2_val in chi2_values]
fig, ax = plt.subplots(figsize=(12, 8))
contour_filled = ax.contourf(N_grid, BR_grid, Z, levels=40, cmap='viridis')
plt.colorbar(contour_filled, ax=ax, label='Joint PDF')
contour_lines = plt.contour(N_grid, BR_grid, Z, levels=sorted(contour_levels), colors=['black', 'black', 'black'])
ax.clabel(contour_lines,  fmt={contour_levels[0]: probability_levels[0], contour_levels[1]: probability_levels[1], contour_levels[2]: probability_levels[2]}, inline=True, fontsize=10)


ax.set_xlabel('Effective value of mosquito population')
ax.set_ylabel('Effective value of basic reproductive number')
ax.set_title('Joint pdf of initiall spillover, Orange County, California')
ax.set_yticklabels([])
plt.grid(True)
plt.xlim(0, 1.3 * max(N_array))
plt.ylim(0, 1.3 * max(BR_array))
plt.savefig('initiall-jointpdforange.pdf', format='pdf', dpi=1000)
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
plt.savefig('JointinitialPDForange', format='pdf', dpi=1000)
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
#years = [year for year in range(1981, 2024)]
years = [year for year in range(2006, 2024) if year!=2010]


ttt=0

for year in years:
    if year>=1981 and year<=2000:
        Casetime=list(np.loadtxt(f"./orangecases{2006}.txt"))
        C=list(np.loadtxt(f"./Cf{year}.txt"))
        d0=get_temperature_for_year(file_path, year)[0]
        dama=[d0[i] for i in range(0,365)]
    elif year>=2000 and year<=2005:
        Casetime=list(np.loadtxt(f"./orangecases{2006}.txt"))
        C=list(np.loadtxt(f"./Cf{year}.txt"))
        dama = list(np.loadtxt(f"./T{year}.txt"))
    elif year==2021:
        Casetime=list(np.loadtxt(f"./orangecases{year}.txt"))
        C=list(np.loadtxt(f"./Cf2013.txt"))
        # C=list(np.loadtxt("./C_esti_2021.txt"))
        dama = list(np.loadtxt(f"./T{year}.txt"))
    elif year==2022:
        Casetime=list(np.loadtxt(f"./orangecases{year}.txt"))
        C=list(np.loadtxt(f"./Cf2013.txt"))
        # C=list(np.loadtxt("./C_esti_2022.txt"))
        dama = list(np.loadtxt(f"./T{year}.txt"))
    elif year==2023:
        Casetime=list(np.loadtxt(f"./orangecases{year}.txt"))
        C=list(np.loadtxt(f"./Cf2013.txt"))
        # C=list(np.loadtxt("./C_esti_2023.txt"))
        dama = list(np.loadtxt(f"./T{year}.txt"))
    else:
        Casetime=list(np.loadtxt(f"./orangecases{year}.txt"))
        dama = list(np.loadtxt(f"./T{year}.txt"))
        C=list(np.loadtxt(f"./Cf{year}.txt"))
    
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
    t=np.linspace(0,365,364)
    x0=(3000000,0,0,0,10,10,10,10,10,100,100,100,100,100,100)
    x=odeint(BMPPT,x0,t)
    muB=(0.15/365)
    incB=1/3
    rC=1/6
    muwnvB=0.9
    R=x[0:365,3]
    I=x[0:365,2]
    M=np.cumsum(x[0:365,6])
    Mos=x[0:365,6]
    BS=x[0:365,11]
    TH=x[0:365,0]+x[0:365,1]+x[0:365,2]+x[0:365,3]
    D=x[0:365,11]+x[0:365,12]+x[0:365,13]+x[0:365,14]
    BetaBM=Coematrix[7,0:364]*Coematrix[3,0:364]/(D+TH)
    BetaMB=Coematrix[3,0:364]/(D+TH)
    R2=(BetaBM*BetaMB*Mos*BS*Coematrix[6,0:364]*incB)/((muB+incB)*(rC+muwnvB+muB)*(Coematrix[2,0:364])*(Coematrix[6,0:364]+Coematrix[2,0:364]))
    R0=[np.sqrt(R2[i]) for i in range(0,364)]
    basicR0=R0
    R0=np.cumsum(R0)
    
    
#----------------------------------------savinf result as txt (mosquitoe profiles)
    # df = pd.DataFrame(M, columns=['Numbers'])
    # df.to_csv(f"./Mosquitoe_Profile_Orange_County_{year}.txt", index=False, header=False, sep='\t')
#------------------------------------------------------------------------------
#--------------------------------------------------------------------------------------
    fig, ax1 = plt.subplots(figsize=(10, 6))
        # Helper function
    def is_inside(collection, x, y):
        """Return True if (x,y) is inside any of the paths in the given contour collection."""
        for path in collection.get_paths():
            if path.contains_point((x, y)):
                return True
        return False
    
    # Identify which collection corresponds to which ellipse
    col_largest  = contour_lines.collections[0]  # 95% ellipse
    col_middle   = contour_lines.collections[1]  # 97% ellipse
    col_smallest = contour_lines.collections[2]  # 99% ellipse
    
    for i in range(364):
        x, y = M[i], R0[i]
        
        # Determine color based on the ellipse tests
        if is_inside(col_smallest, x, y):
            color = 'r'
        elif is_inside(col_middle, x, y):
            color = 'orange'
        elif is_inside(col_largest, x, y):
            color = 'y'
        else:
            color = 'g'
    
        # Old: a narrow line at x=i
        # plt.axvline(x=i, color=color, linestyle='solid')
        
        # New: fill from day i to day i+1
        plt.axvspan(i, i + 1, facecolor=color, edgecolor=color)

    
       
        plt.axvline(x=i, color=color, linestyle='solid')
    if year > 2020:
        
    
        # Plot mosquito profile on the left y-axis
        ax1.plot(Mos / max(Mos), '-b', label='Predicted profile of mosquitoes',linewidth=3)
        ax1.set_xlabel('Days')
        ax1.set_ylabel('Predicted relative abundance of mosquitoes', color='b')
        ax1.tick_params(axis='y', labelcolor='b')
        ax1.set_title('Predicted risk - Orange County, California - ' + str(year))
    
        # Plot spillover time
        
    
        # Create second y-axis for basic reproductive number
        # ax2 = ax1.twinx()
        # ax2.plot(moving_average(basicR0, window_size=8), '-k', label='Predicted basic reproductive number',linewidth=3)
        ax1.axvline(s, color='k', linestyle='--', label='First spillover time')
        # ax2.set_ylabel('Predicted basic reproductive number', color='k')
        # ax2.tick_params(axis='y', labelcolor='k')
    
        # Add legends for both plots
        fig.legend(loc='upper left', bbox_to_anchor=(0.1, 0.8))
        fig.tight_layout()
    
        # Save the figure
        plt.savefig(f'predicted_initial_risk_{year}.pdf', format='pdf', dpi=1000)
        plt.show()
    
    else:

    
        # Plot mosquito profile on the left y-axis
        ax1.plot(Mos / max(Mos), '-b', label='Estimated profile of mosquitoes',linewidth=3)
        ax1.set_xlabel('Days')
        ax1.set_ylabel('Estimated relative abundance of mosquitoes', color='b')
        ax1.axvline(s, color='k', linestyle='--', label='First spillover time')
        ax1.tick_params(axis='y', labelcolor='b')
        
        ax1.set_title('Estimated risk - Orange County, California - ' + str(year))
    
        # Plot spillover time
        
    
        # Create second y-axis for basic reproductive number
        # ax2 = ax1.twinx()
        # ax2.plot(moving_average(basicR0, window_size=8), '-k', label='Estimated basic reproductive number',linewidth=3)
        ax1.axvline(
            s, 
            color='k', 
            linestyle='--', 
            linewidth=3,   # <--- increase thickness here
            label='Initial spillover time'
        )        # ax2.set_ylabel('Estimated basic reproductive number', color='k')
        # ax2.tick_params(axis='y', labelcolor='k')
    
        # Add legends for both plots
        if year==2012:
            fig.legend(loc='upper left', bbox_to_anchor=(0.1, 0.8))
        fig.tight_layout()
    
        # Save the figure
        plt.savefig(f'estimated_initial_risk_{year}.pdf', format='pdf', dpi=1000)
    plt.show()
    
    
    