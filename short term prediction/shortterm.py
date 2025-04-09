# -*- coding: utf-8 -*-
"""
Created on Thu Apr 25 12:58:31 2024

@author: shosseini
"""


import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from statsmodels.tsa.ar_model import AutoReg
from scipy.integrate import odeint
from scipy.stats import gaussian_kde




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
T2021 = list(np.loadtxt("./T2021.txt"))
H2021 = list(np.loadtxt("./H2021.txt"))
Cf2021 = list(np.loadtxt("./Cf2021.txt"))
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




years = [2006, 2007, 2008, 2009, 2011, 2012, 2013, 2014, 2015, 2016, 2017, 2018, 2019, 2020,2021]

# Load all C{year} files into a 2D NumPy array
C_values = np.array([np.loadtxt(f"./Cf{year}.txt") for year in years])

# Compute the element-wise mean across all years
meanC = np.mean(C_values, axis=0)

# Convert to a list if needed
meanC = meanC.tolist()

print("mean",meanC)







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
targetyear=2023




















years = [year for year in range(2006, 2022) if year != 2010]
for year in years:
    redyear=[]
    oranyear=[]
    yellyear=[]
    redtogreen=[]
    greenyear=[]
    risktoyear=[]
    relativeredtorisk=[]
    dama = list(np.loadtxt(f"./T{year}.txt"))
    if year <2022:
        C=list(np.loadtxt(f"./Cf{year}.txt"))
    else:
        C=meanC#list(np.loadtxt(f"./C_esti_{year}.txt"))
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
                return [i] #* int(value)  # Convert to int
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
    print("ct",ct,year)
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

   
    for i in ct:
        #print(M[i*7])
        B.append(R0[i*7])
        N.append(M[i*7])


# ------------------------------------------------------------------
# Assuming you already have B and N as NumPy arrays
# B = np.array(B)
# N = np.array(N)

# 1) Combine B and N into a 2D array (shape: (2, n_samples))
data = np.vstack([N, B])

# 2) Fit a 2D Gaussian KDE
kde = gaussian_kde(data, bw_method=0.2)

# ------------------------------------------------------------------
# 3) Define the domain for N and B, and generate a grid
n_min, n_max = 1e12, 5e12
b_min, b_max = 70, 220

n_points = 100  # number of points in each direction
n_range = np.linspace(n_min, n_max, n_points)
b_range = np.linspace(b_min, b_max, n_points)
N_grid, B_grid = np.meshgrid(n_range, b_range)

# Evaluate the KDE on the grid
positions = np.vstack([N_grid.ravel(), B_grid.ravel()])
Z = kde(positions).reshape(N_grid.shape)

# ------------------------------------------------------------------
# 4) Normalize the density over the grid to ensure we can treat it as a discrete PDF
#
#    The sum of Z over the grid times the area of each cell (dx*dy) should be ~1.
#    We'll use this approach to find the fraction of total probability.

dx = (n_max - n_min) / (n_points - 1)  # width of each cell in N direction
dy = (b_max - b_min) / (n_points - 1)  # height of each cell in B direction
cell_area = dx * dy

# Flatten Z for sorting in descending order
Z_flat = Z.ravel()
# Sort densities high->low
sorted_indices = np.argsort(Z_flat)[::-1]
Z_sorted = Z_flat[sorted_indices]

# Cumulative sum of sorted densities times the area of each cell
cdf = np.cumsum(Z_sorted) * cell_area
# Normalize so that total probability ~ 1
cdf /= cdf[-1]

# ------------------------------------------------------------------
# 5) Find the density thresholds for 90%, 95%, and 99%
def find_level(prob):
    """
    Returns the density threshold z* at which 'prob' fraction of total mass is enclosed.
    """
    # Find the index where cdf crosses 'prob'
    idx = np.searchsorted(cdf, prob)
    return Z_sorted[idx]
a=0.8
b=0.85
c=0.9
# Compute probability levels
level_90 = find_level(a)
level_95 = find_level(b)
level_99 = find_level(c)

# ------------------------------------------------------------------
# 6) Plot the 2D PDF using contourf
fig, ax = plt.subplots(figsize=(8, 6))

# Filled contours for the underlying density
contourf_plot = ax.contourf(N_grid, B_grid, Z, levels=50, cmap='viridis')
cbar = plt.colorbar(contourf_plot, ax=ax)
cbar.set_label('PDF (Gaussian KDE)')

# ------------------------------------------------------------------
# 7) Overlay the contour lines at the 90%, 95%, 99% thresholds
levels = [level_90, level_95, level_99]
contours = ax.contour(N_grid, B_grid, Z, levels=sorted(levels), colors=['black', 'red', 'blue'])

# Label the contours with their respective probability
fmt = {
    min(levels): '70%',
    np.median(levels): '75%',
    max(levels): '85%'
}
ax.clabel(contours, inline=True, fontsize=10, fmt=fmt)

# ------------------------------------------------------------------
# 8) Modify axis settings
ax.set_xlabel('M')  # X-axis label remains
ax.set_ylabel('R')   # Remove Y-axis label
ax.set_yticklabels([])  # Remove Y-axis tick labels

# Add title
ax.set_title('2D PDF of (M, R)')

# Remove extra padding and save as PDF
plt.tight_layout()
#plt.savefig("initiall-jointpdforange.pdf", format="pdf", dpi=1000, bbox_inches="tight")

# Show plot (optional)
plt.show()






# -------------------------------------------------------------------------
# 1. Fit your 2D KDE to arrays B and N.

# Suppose B and N are already defined as NumPy arrays:
#   B = np.array([...])
#   N = np.array([...])

# Combine into shape (2, n_samples) for gaussian_kde
data = np.vstack([N, B])

# Fit the KDE
kde = gaussian_kde(data, bw_method=0.2)

# -------------------------------------------------------------------------
# 2. Build a grid covering the domain, evaluate the KDE, and compute 
#    the cumulative distribution of density values in descending order.

# Define a reasonable bounding box for your data
n_min, n_max = 0.5e12, 5e12
b_min, b_max = 50, 250

n_points = 100  # increase for better accuracy
n_range = np.linspace(n_min, n_max, n_points)
b_range = np.linspace(b_min, b_max, n_points)
N_grid, B_grid = np.meshgrid(n_range, b_range)

# Evaluate KDE on this grid
positions = np.vstack([N_grid.ravel(), B_grid.ravel()])
Z = kde(positions).reshape(N_grid.shape)

# Calculate area of each cell in the grid (for discrete integration)
dx = (n_max - n_min) / (n_points - 1)
dy = (b_max - b_min) / (n_points - 1)
cell_area = dx * dy

# Flatten the density array, sort in descending order
Z_flat = Z.ravel()
sort_idx = np.argsort(Z_flat)[::-1]     # indices of Z_flat, largest -> smallest
Z_sorted = Z_flat[sort_idx]

# Compute the cumulative sum (descending density)
cdf_desc = np.cumsum(Z_sorted) * cell_area
# Normalize so total probability is ~ 1
cdf_desc /= cdf_desc[-1]

# -------------------------------------------------------------------------
# 3. Define a helper function that, given a density d, finds
#    p(d) = Prob[f(X) >= d].

def fraction_enclosed(d):

    idx = np.searchsorted(Z_sorted[::-1], d, side='right')
    
    i_descending = len(Z_sorted) - idx  # this is the index in Z_sorted (descending)
    
    if i_descending <= 0:
        # That means d is bigger than any density in Z_sorted
        # => fraction is 0 (no region has density >= d).
        return 0.0
    elif i_descending >= len(Z_sorted):
        # d is smaller than or equal to the smallest density
        # => entire domain has density >= d => fraction is 1.
        return 1.0
    else:
        return cdf_desc[i_descending - 1]

# -------------------------------------------------------------------------
# 4. Define the classification function using fraction_enclosed(density).
#    We want four disjoint shells: inside 90%, between 90–95%, 
#    between 95–99%, or outside 99%.

def is_inside(n_value, b_value):
    """
    Classify the point (n_value, b_value) into one of four disjoint shells:
      - "Inside 90%"         
      - "Between 90% and 95%"
      - "Between 95% and 99%"
      - "Outside 99%"
    based on the fraction of the KDE's probability mass that lies 
    at densities >= f(n_value,b_value).
    """
    # Evaluate the KDE at (n_value, b_value)
    d = kde(np.array([[n_value], [b_value]]))[0]
    
    # Compute p(d) = Probability of having density >= d
    p = fraction_enclosed(d)
    
    # Now classify according to p
    if p <= a:
        return 0#"Inside 90%"
    elif p <= b:
        return 1#"Between 90% and 95%"
    elif p <=c:
        return 2#"Between 95% and 99%"
    else:
        return 3#"Outside 99%"

# -------------------------------------------------------------------------
# 5. Example Usage

test_n = 1.8e12   # some test point in your domain
test_b = 120
region = is_inside(test_n, test_b)
print(f"The point (N={test_n}, B={test_b}) is classified as: {region}")





#------------------------------------------------------------------------------
#-------------------------------------prediction part--------------------------
# thid code will return the risk prediction for two weeks 
# the vector of the carying capacity (as well as temperature) should be importing to
# the code (so using another code the carying capacity should be predicted.)
# at the end it is plotted based in the perioi of iteration you want.



# Slice the dama list to obtain the temperature data for each period
if targetyear == 2022:
    dama = list(np.loadtxt('2000-2022tem.txt'))
elif targetyear ==2023:
    dama = list(np.loadtxt('2000-2023tem.txt'))
elif targetyear == 2024:
    dama = list(np.loadtxt('2000-2024tem.txt'))












































T=dama
prediction_lead=14
dama=dama[len(dama)-2*365:len(dama)]
for datenumber in range(150,365,5):#365,1):
#for datenumber in {345}:
    # datenumber=260
    d=365-datenumber
    inuse=dama[0:len(dama)-d]
    #print(inuse)
    model = AutoReg(inuse, lags=30)
    model_fitted = model.fit()
    predictions = model_fitted.predict(start=len(inuse), end=len(inuse)+prediction_lead)
    T2=list(predictions)
       
    temp=T[len(T)-365:len(T)-365+datenumber]
    for i in range(0,len(predictions)):
        temp.append(np.round(predictions[i],decimals=2))

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
    
    
    C =list(np.loadtxt(f"./C_esti_{targetyear}.txt"))
    Casetime=list(np.loadtxt(f"./orangecases{targetyear}.txt"))
        
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
    muB=(0.15/365)
    incB=1/3
    rC=1/6
    muwnvB=0.9
    R=x[0:datenumber+1+prediction_lead,3]
    I=x[0:datenumber+1+prediction_lead,2]
    M=x[0:datenumber+1+prediction_lead,6]
    MM=np.cumsum(M)
    Mos=M
    BS=x[0:datenumber+1+prediction_lead,11]
    TH=x[0:datenumber+1+prediction_lead,0]+x[0:datenumber+1+prediction_lead,1]+x[0:datenumber+1+prediction_lead,2]+x[0:datenumber+1+prediction_lead,3]
    D=x[0:datenumber+1+prediction_lead,11]+x[0:datenumber+1+prediction_lead,12]+x[0:datenumber+1+prediction_lead,13]+x[0:datenumber+1+prediction_lead,14]
    BetaBM=Coematrix[7,0:datenumber+prediction_lead]*Coematrix[3,0:datenumber+prediction_lead]/(D+TH)
    BetaMB=Coematrix[3,0:datenumber+prediction_lead]/(D+TH)
    R2=(BetaBM*BetaMB*Mos*BS*Coematrix[6,0:datenumber+prediction_lead]*incB)/((muB+incB)*(rC+muwnvB+muB)*(Coematrix[2,0:datenumber+prediction_lead])*(Coematrix[6,0:datenumber+prediction_lead]+Coematrix[2,0:datenumber+prediction_lead]))
    R0=[np.sqrt(R2[i]) for i in range(0,datenumber+prediction_lead)]
    basicR0=R0
    R0=np.cumsum(R0)
    

    
    
    RR=0
    Or=[]
    Yel=[]
    RED=[]
    plt.figure(figsize=(8, 5))
    for i in range(0,datenumber+prediction_lead):
        x, y = MM[i], R0[i]
        
        # Determine color based on the ellipse tests
        if is_inside( x, y)==0:
            color = 'r'
            #plt.axvline(i, color='red')
            RED.append(1)
        elif is_inside(x, y)==1:
            color = 'orange'
            #plt.axvline(i, color='orange')
            Or.append(1)
        elif is_inside(x, y)==2:
            color = 'y'
           # plt.axvline(i, color='y')
            Yel.append(1)
        else:
            color = 'g'
            plt.axvline(i, color='g')
            x, y = M[i], R0[i]
        
            
        plt.axvspan(i, i + 1, facecolor=color, edgecolor=color)

    

    # plt.title('Probability of spillover for 200'+str(2021))
    # plt.xlabel('day')
    # plt.legend(['green:P < 0.03', 'yellow: 0.03 < P < 0.3', 'orange: 0.3 < P < 0.6', 'red: 0.6 < P '],
    #           loc='upper left')
    #plt.savefig('200' + str(f'{r}') + '.pdf')
    # plt.savefig('200' + str(f'{year}'))
    plt.ylim(0, 12**10)
    plt.xlim(0,370)
    plt.yticks([])
    plt.plot(M,'--w',label='Estimated profile of mosquitoes')
    #plt.axvline(s, color='b', linestyle='--', label='spillover time')

    #print(ct*7)
    for i in ct:
        if i==min(ct):
            plt.axvline(x=i*7, color='b', linestyle='--', label='Spillover time')
    plt.title(f'Two weeks prediction,Orange County {targetyear}', fontsize=15)
    plt.xlabel('Day', fontsize=14)
    plt.axvline(x=datenumber, color='k', linestyle='--', linewidth=1, label='Prediction period')
    plt.axvline(x=datenumber+prediction_lead, color='k', linestyle='--', linewidth=1)
    if datenumber==195:
        plt.legend(loc='upper left', fontsize='small')  # You can also use specific sizes like 10, 12, etc.
    plt.savefig(f'short term-oregan-risk{datenumber}.pdf', format='pdf', dpi=1000)
    plt.show()
    # if RR>0:
    #     print( '# red day', RR, '# orange day', oo, '# of yellow day', yy)
print(sum(RED),sum(RED)+sum(Or)+sum(Yel))
    
  
