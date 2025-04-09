# -*- coding: utf-8 -*-
"""
Created on Mon Feb 17 23:56:12 2025

@author: shosseini
"""

import numpy as np
import matplotlib.pyplot as plt
from scipy.integrate import odeint
from scipy.stats import gamma
from scipy.stats import gaussian_kde
import pandas as pd
from scipy.optimize import minimize
from scipy.stats import gamma, norm, multivariate_normal, chi2, kstest
from collections import Counter
from scipy.interpolate import RegularGridInterpolator
from matplotlib.colors import Normalize
from matplotlib.cm import ScalarMappable
from shapely.geometry import LineString, Point
import os
from statsmodels.tsa.ar_model import AutoReg
from scipy.stats import poisson
#--------------------------------------Section of importing data 
C2006 = list(np.loadtxt("./Cf2006.txt"))
#T2006 = list(np.loadtxt("./T2006.txt"))
# H2006 = list(np.loadtxt("./H2006.txt"))
#---------------------------------------------------------------
C2007 = list(np.loadtxt("./Cf2007.txt"))
# #T2007 = list(np.loadtxt("./T2007.txt"))
# H2007 = list(np.loadtxt("./H2007.txt"))
#---------------------------------------------------------------
C2008 = list(np.loadtxt("./Cf2008.txt"))
# T2008 = list(np.loadtxt("./T2008.txt"))
# H2008 = list(np.loadtxt("./H2008.txt"))
#---------------------------------------------------------------
C2009 = list(np.loadtxt("./Cf2009.txt"))
# T2009 = list(np.loadtxt("./T2009.txt"))
# H2009 = list(np.loadtxt("./H2009.txt"))
#---------------------------------------------------------------
C2011 = list(np.loadtxt("./Cf2011.txt"))
# T2011 = list(np.loadtxt("./T2011.txt"))
# H2011 = list(np.loadtxt("./H2011.txt"))
#---------------------------------------------------------------
C2012 = list(np.loadtxt("./Cf2012.txt"))
# T2012 = list(np.loadtxt("./T2012.txt"))
# H2012 = list(np.loadtxt("./H2012.txt"))
#---------------------------------------------------------------
C2013 = list(np.loadtxt("./Cf2013.txt"))
# T2013 = list(np.loadtxt("./T2013.txt"))
# H2013 = list(np.loadtxt("./H2013.txt"))
#---------------------------------------------------------------
C2014 = list(np.loadtxt("./Cf2014.txt"))
# T2014 = list(np.loadtxt("./T2014.txt"))
# H2014 = list(np.loadtxt("./H2014.txt"))
#---------------------------------------------------------------
C2015 = list(np.loadtxt("./Cf2015.txt"))
# T2015 = list(np.loadtxt("./T2015.txt"))
# H2015 = list(np.loadtxt("./H2015.txt"))
#---------------------------------------------------------------
C2016 = list(np.loadtxt("./Cf2016.txt"))
# T2016 = list(np.loadtxt("./T2016.txt"))
# H2016 = list(np.loadtxt("./H2016.txt"))
#---------------------------------------------------------------
C2017 = list(np.loadtxt("./Cf2017.txt"))
# T2017 = list(np.loadtxt("./T2017.txt"))
# H2017 = list(np.loadtxt("./H2017.txt"))
#---------------------------------------------------------------
C2018 = list(np.loadtxt("./Cf2018.txt"))
# T2018 = list(np.loadtxt("./T2018.txt"))
# H2018 = list(np.loadtxt("./H2018.txt"))
#---------------------------------------------------------------
C2019 = list(np.loadtxt("./Cf2019.txt"))
# T2019 = list(np.loadtxt("./T2019.txt"))
# H2019 = list(np.loadtxt("./H2019.txt"))
#---------------------------------------------------------------
C2020 = list(np.loadtxt("./Cf2020.txt"))
# T2020 = list(np.loadtxt("./T2020.txt"))
# H2020 = list(np.loadtxt("./H2020.txt"))
#---------------------------------------------------------------


#====================================functions=============================

def get_temperature_for_year(file_path, year):
    """
    Read the CSV file and get the temperature data for the specified year.
    
    Parameters:
        file_path (str): Path to the CSV file.
        year (int): Year for which to extract temperature data.
    
    Returns:
        list: List of computed temperature values for the given year.
    """
    # Read the CSV file
    df = pd.read_csv(file_path)
    
    # Filter the data for the specified year
    df_year = df[df['Year'] == year]
    
    # Convert the columns to lists
    p = df_year['P'].tolist()
    t = df_year['T'].tolist()
    h = df_year['H'].tolist()
    
    # Compute the element-wise average
    result = [ 0.5*(pi + ti + hi) for pi, ti, hi in zip(p, t, h)]
    
    return result,t,h,p


# File path and years to process
file_path = 'infoupdated.csv'
#Latitude and longitude coordinates are: 33.787914, -117.853104

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



#===============================dictionary for saving number of cases and returning year and week of that
#=========================================================================================================
#=========================================================================================================

years = list(range(2006, 2022))  # Years from 2006 to 2021 (excluding 2022)

# Dictionary to store weeks and years for each number of cases
cases_occurrences = {}

# Read case data for each year
for year in years:
    try:
        # Load case data for the year
        case_data = np.loadtxt(f"./orangecases{year}.txt")
        
        # Loop through each week and count cases
        for week, cases in enumerate(case_data, start=1):  # Week is 1-based
            if cases > 0:  # Only consider weeks with reported cases
                cases = int(cases)  # Convert cases to integer
                if cases not in cases_occurrences:
                    cases_occurrences[cases] = []  # Initialize list for this case count
                cases_occurrences[cases].append((week-1, year))  # Add (week, year) pair

    except FileNotFoundError:
        print(f"Data for {year} not found. Skipping this year.")


#=========================================================================================================
#=========================================================================================================
#================================================Mean of carrying capacities=============================

years = [2006, 2007, 2008, 2009, 2011, 2012, 2013, 2014, 2015, 2016, 2017, 2018, 2019, 2020]

# Load all C{year} files into a 2D NumPy array
C_values = np.array([np.loadtxt(f"./Cf{year}.txt") for year in years])

# Compute the element-wise mean across all years
meanC = np.mean(C_values, axis=0)

# Convert to a list if needed
meanC = meanC.tolist()

#------------------------------------------------------------------------------
y=2022
file_path = 'infoupdated.csv'
dama=get_temperature_for_year(file_path, y)[1]
C=meanC



# Placeholder for coefficient matrix
Coematrix = np.zeros((8, 365))  

# Fill Coematrix based on time-dependent functions

    
for time in range(0, 365):
    Coematrix[0, time] = gammaMP(time)
    Coematrix[1, time] = pLSP(time)
    Coematrix[2, time] = muMP(time)
    Coematrix[3, time] = BRP(time)
    Coematrix[4, time] = alphaEP(time)
    Coematrix[5, time] = bEP(time)
    Coematrix[6, time] = PDRP(time)
    Coematrix[7, time] = VCP(time)

# Define BMPPT function
def BMPPT(x, time):
    time = int(divmod(time, 365)[1])
    x = np.maximum(x, 0)  # Ensure non-negative states

    # Constants and parameters
    FB = 2
    EVB = 0.021
    muE = 4.9315 * (10**-4)
    BDR = 9.041 * (10**-4)
    muF = 6.30136 * (10**-4)
    muB = (0.15 / 365)
    incB = 1 / 3
    recoveryrateC = 1 / 6
    muwnvB = 0.9
    muHbirth = 5
    muhumanWND = 0.01
    incHumanWND = 1 / 6
    infectionhumanWND = 1 / 3
    muHuman = (1 / (77 * 365))

    TH = x[0] + x[1] + x[2] + x[3]
    D = x[11] + x[12] + x[13] + x[14]
    ph = Coematrix[3, time] / (D + TH)

    # Human equations
    dsHdt = (muHuman * (TH)) - ph * Coematrix[7, time] * x[8] * x[0] - muHuman * x[0]
    deHdt = ph * Coematrix[7, time] * x[8] * x[0] - (incHumanWND) * x[1] - muHuman * x[1]
    diHdt = (incHumanWND) * x[1] - (infectionhumanWND) * x[2] - (muhumanWND) * x[2] - muHuman * x[2]
    drHdt = (infectionhumanWND) * x[2] - muHuman * x[3]

    # Mosquito equations
    deggMdt = Coematrix[4, time] * (x[6] + x[7] + x[8]) - Coematrix[5, time] * x[4]
    dacuMdt = Coematrix[5, time] * x[4] * max(0, (1 - (x[5] / (C[time])))) - Coematrix[0, time] * Coematrix[1, time] * x[5]
    dsMdt = Coematrix[0, time] * Coematrix[1, time] * x[5] - ph * x[13] * x[6] - Coematrix[2, time] * x[6]
    deMdt = ph * x[13] * x[6] - Coematrix[6, time] * x[7] - Coematrix[2, time] * x[7]
    diMdt = Coematrix[6, time] * x[7] - Coematrix[2, time] * x[8]

    # Bird equations
    deggBdt = FB * (x[11] + x[12] + x[13] + x[14]) - EVB * x[9] - muE * x[9]
    dfleBdt = EVB * x[9] * max(0, (1 - (x[10] / 20000))) - BDR * x[10] - muF * x[10]
    dsBdt = BDR * x[10] - ph * Coematrix[7, time] * x[8] * x[11] - muB * x[11]
    deBdt = ph * Coematrix[7, time] * x[8] * x[11] - incB * x[12] - muB * x[12]
    diBdt = incB * x[12] - muwnvB * x[13] - muB * x[13] - recoveryrateC * x[13]
    drBdt = recoveryrateC * x[13] - muB * x[14]

    return (
        dsHdt, deHdt, diHdt, drHdt,
        deggMdt, dacuMdt, dsMdt, deMdt, diMdt,
        deggBdt, dfleBdt, dsBdt, deBdt, diBdt, drBdt
    )

    
 
t = np.linspace(0, 365, 364)
x0 = (3000000, 0, 0, 0, 10, 10, 10, 10, 10, 100, 100, 100, 100, 100, 100)
x = odeint(BMPPT, x0, t)
incB=1/3
rC=1/6
muwnvB=0.9
muB=(0.15/365)

# Compute cumulative mosquito values (M)

M=np.cumsum(x[0:365,6])
Mos=x[0:365,6]




# Assume Mos is your array
np.savetxt("./MMosquito_Profile_2022.txt", Mos, fmt="%.6f")  # Save with 6 decimal precision


BR_dic1=np.cumsum((get_temperature_for_year(file_path, y)[0]))

Curve = []
for i in range(0, 364):
    #print(i)
    Curve.append((M[i], BR_dic1[i]))
#print(Curve)
 
# Save to a text file in two columns format
with open("Curve-2022.txt", "w") as file:
    for item in Curve:
        file.write(f"{item[0]}\t{item[1]}\n")













# List of years to process








Time=list(np.loadtxt('./orangecases2022.txt'))    # vector of cases to find out atarting point
plt.plot(Time)
plt.show()

def first_nonzero_index(lst):
    """
    Returns the index of the first nonzero element in lst,
    or None if no such element exists.
    """
    for i, val in enumerate(lst):
        if val != 0:
            return i
    return None


starting_time=int(first_nonzero_index(Time)*7)
starting=starting_time



tt=0
accuracy=0.1









years = [year for year in range(2006, 2022) if year != 2010]
print("year",years)
# Dictionary to store cumulative mosquito profiles (M) and weather parameter (W) by year
mosquito_profiles = {}
BR_dic={}


#d=365
for year in years:
    if year==2021:
        C=meanC
    else:
        C = list(np.loadtxt(f"./Cf{year}.txt"))
    dama = get_temperature_for_year(file_path, year)[1]
    Casetime = list(np.loadtxt(f"./orangecases{year}.txt"))

    # Placeholder for coefficient matrix
    Coematrix = np.zeros((8, 365))  

    # Fill Coematrix based on time-dependent functions

        
    for time in range(0, 365):
        Coematrix[0, time] = gammaMP(time)
        Coematrix[1, time] = pLSP(time)
        Coematrix[2, time] = muMP(time)
        Coematrix[3, time] = BRP(time)
        Coematrix[4, time] = alphaEP(time)
        Coematrix[5, time] = bEP(time)
        Coematrix[6, time] = PDRP(time)
        Coematrix[7, time] = VCP(time)

    # Define BMPPT function
    def BMPPT(x, time):
        time = int(divmod(time, 365)[1])
        x = np.maximum(x, 0)  # Ensure non-negative states

        # Constants and parameters
        FB = 2
        EVB = 0.021
        muE = 4.9315 * (10**-4)
        BDR = 9.041 * (10**-4)
        muF = 6.30136 * (10**-4)
        muB = (0.15 / 365)
        incB = 1 / 3
        recoveryrateC = 1 / 6
        muwnvB = 0.9
        muHbirth = 5
        muhumanWND = 0.01
        incHumanWND = 1 / 6
        infectionhumanWND = 1 / 3
        muHuman = (1 / (77 * 365))

        TH = x[0] + x[1] + x[2] + x[3]
        D = x[11] + x[12] + x[13] + x[14]
        ph = Coematrix[3, time] / (D + TH)

        # Human equations
        dsHdt = (muHuman * (TH)) - ph * Coematrix[7, time] * x[8] * x[0] - muHuman * x[0]
        deHdt = ph * Coematrix[7, time] * x[8] * x[0] - (incHumanWND) * x[1] - muHuman * x[1]
        diHdt = (incHumanWND) * x[1] - (infectionhumanWND) * x[2] - (muhumanWND) * x[2] - muHuman * x[2]
        drHdt = (infectionhumanWND) * x[2] - muHuman * x[3]

        # Mosquito equations
        deggMdt = Coematrix[4, time] * (x[6] + x[7] + x[8]) - Coematrix[5, time] * x[4]
        dacuMdt = Coematrix[5, time] * x[4] * max(0, (1 - (x[5] / (C[time])))) - Coematrix[0, time] * Coematrix[1, time] * x[5]
        dsMdt = Coematrix[0, time] * Coematrix[1, time] * x[5] - ph * x[13] * x[6] - Coematrix[2, time] * x[6]
        deMdt = ph * x[13] * x[6] - Coematrix[6, time] * x[7] - Coematrix[2, time] * x[7]
        diMdt = Coematrix[6, time] * x[7] - Coematrix[2, time] * x[8]

        # Bird equations
        deggBdt = FB * (x[11] + x[12] + x[13] + x[14]) - EVB * x[9] - muE * x[9]
        dfleBdt = EVB * x[9] * max(0, (1 - (x[10] / 20000))) - BDR * x[10] - muF * x[10]
        dsBdt = BDR * x[10] - ph * Coematrix[7, time] * x[8] * x[11] - muB * x[11]
        deBdt = ph * Coematrix[7, time] * x[8] * x[11] - incB * x[12] - muB * x[12]
        diBdt = incB * x[12] - muwnvB * x[13] - muB * x[13] - recoveryrateC * x[13]
        drBdt = recoveryrateC * x[13] - muB * x[14]

        return (
            dsHdt, deHdt, diHdt, drHdt,
            deggMdt, dacuMdt, dsMdt, deMdt, diMdt,
            deggBdt, dfleBdt, dsBdt, deBdt, diBdt, drBdt
        )

        
     
    t = np.linspace(0, 365, 364)
    x0 = (3000000, 0, 0, 0, 10, 10, 10, 10, 10, 100, 100, 100, 100, 100, 100)
    x = odeint(BMPPT, x0, t)
    incB=1/3
    rC=1/6
    muwnvB=0.9
    muB=(0.15/365)

    # Compute cumulative mosquito values (M)
    M=np.cumsum(x[0:365,6])
    Mos=x[0:365,6]

    
    mosquito_profiles[year] =M
    BR_dic[year]=np.cumsum((get_temperature_for_year(file_path, year)[0]))



data=[]
dataBR=[]

data1=[]
dataBR1=[]

data2=[]
dataBR2=[]



border=13
A=[]
B=[]
for i in range(0,32):
    num_cases = i
    #print(f"Processing for {num_cases} cases")
    
    # Check if the num_cases key exists in the dictionary
    if num_cases in cases_occurrences:
        occurrences = cases_occurrences[num_cases]
        
        # Dictionary to save cumulative mosquito values for each occurrence
        mosquito_values_list = []
        Bird_value_list=[]
        
        # Process the occurrences and save mosquito values
        for week, year in occurrences:
            if year in mosquito_profiles:
                # for k in range(0,1):                    
                #     mosquito_value = mosquito_profiles[year][7 * (week-1-k)]  # Week is 1-based, adjust for 0-indexed array
                #     BR_value=BR_dic[year][7 * (week-1-k)]
                if i<border:
                    for k in range(0,1):
                        mosquito_value = mosquito_profiles[year][7 * (week-1-k)]  # Week is 1-based, adjust for 0-indexed array
                        BR_value=BR_dic[year][7 * (week-1-k)]
                        data1.extend([mosquito_value]*num_cases)
                        dataBR1.extend([BR_value]*num_cases)
                    
                else:
                    for k in range(0,1):
                        mosquito_value = mosquito_profiles[year][7 * (week-1-k)]  # Week is 1-based, adjust for 0-indexed array
                        BR_value=BR_dic[year][7 * (week-1-k)]
                        data2.extend([mosquito_value]*num_cases)
                        dataBR2.extend([BR_value]*num_cases)
                        A.append(mosquito_value)
                        B.append(BR_value)
            mosquito_values_list.append(mosquito_value)
            Bird_value_list.append(BR_value)



#--------------------------------------------------------------------------------------

# Combine datasets into coordinate pairs
coordinates1 = list(zip(data1, dataBR1))
coordinates2 = list(zip(data2, dataBR2))

# Combine both datasets into one
all_coordinates = coordinates1 + coordinates2

x1, y1 = zip(*all_coordinates)







def kernel_frequency_surface(coordinates, bw=0.3):
    # Count the frequency of each unique coordinate
    coordinate_counts = Counter(coordinates)
    unique_coords = np.array(list(coordinate_counts.keys()))
    weights = np.array(list(coordinate_counts.values()))  # Frequencies as weights

    # Separate the unique coordinates into x and y values
    x_vals = unique_coords[:, 0]
    y_vals = unique_coords[:, 1]

    # Create the data points for kernel smoothing
    points = np.vstack([x_vals, y_vals])
    kde = gaussian_kde(points, weights=weights, bw_method=bw)

    # ------------------------------------------
    # Force the domain to start at 0
    # ------------------------------------------
    x_min = 0.0  
    x_max = 9.1 * 10**12
    y_min = 0
    y_max = 6000

    # Create a grid for evaluation
    x_grid = np.linspace(x_min, x_max, 400)
    y_grid = np.linspace(y_min, y_max, 400)
    grid_x, grid_y = np.meshgrid(x_grid, y_grid)
    grid_positions = np.vstack([grid_x.ravel(), grid_y.ravel()])

    # Evaluate the smoothed frequency surface
    grid_z = kde(grid_positions).reshape(grid_x.shape)

    # Find the maximum density
    maxden = grid_z.max()
    #print("Maximum density (maxden):", maxden)





    # Scale the output by the desired factor
    scaling_factor = 31 / maxden
    grid_z *= scaling_factor  # Scale the height surface
    np.savez("2d_rate_2022_short.npz", X=grid_x, Y=grid_y, Z=grid_z)
    # print("Scaled 2D rate saved as '2d_rate_2022.npz'")



    # Plot the surface
    plt.figure(figsize=(10, 8))

    # Filled contour
    contour_filled = plt.contourf(grid_x, grid_y, grid_z, levels=50, cmap='viridis')
    plt.colorbar(label='Scaled Frequency')

    # Contour lines (on top of the filled contour)
    contour_lines = plt.contour(grid_x, grid_y, grid_z, levels=10, colors='black', linewidths=1)
    plt.clabel(contour_lines, inline=True, fontsize=8)  # Add labels to the contour lines

    # Scatter plot of original points
    plt.scatter(x_vals, y_vals, color='red', alpha=0.7, s=10, label='observations')

    plt.title('Rate of the Poisson Process')
    plt.xlabel('M')
    plt.ylabel('W')

    # Hide the numeric tick labels on the y-axis (but keep the label itself)
    plt.tick_params(axis='y', labelleft=False)

    plt.legend()
    plt.tight_layout()
    
    plt.xlim(1e12,7e12)
    plt.ylim(2000,5000)
    # plt.savefig('rate_plot.pdf', format="pdf", dpi=300)
    parameter = []
    parameter = np.loadtxt("Curve-2022.txt")

    param_x = [point[0] for point in parameter]
    param_y = [point[1] for point in parameter]
    plt.scatter(param_x, param_y, color='black', alpha=0.5, s=12, label='Parameters')


    plt.show()

    return grid_x, grid_y, grid_z



grid_x, grid_y, grid_z = kernel_frequency_surface(all_coordinates, bw=accuracy)

#=========================================prior==========================================
#=========================================prior==========================================
#=========================================prior==========================================
#=========================================prior==========================================



parameter=np.loadtxt('Curve-2022.txt')
def plot_curve_with_fixed_shadow_and_uniform_pdf(coordinates, r=100, grid_res=100):
    """
    Plots a curve with a "shadow" buffer at distance r,
    sets up a uniform PDF in that region, and returns:
      1) area_of_region (float),
      2) shadow_polygon (Shapely Polygon),
      3) pdf_function (callable: pdf(x,y)).

    Also saves a discrete representation of the PDF in "2d_prior_test.npz":
      - area  (float),
      - X, Y  (2D meshgrid coordinates),
      - Z     (PDF values at each mesh point).

    Args:
        coordinates (list of (float, float)): The (x,y) points of the curve.
        r (float): Radius of the buffer around the curve.
        grid_res (int): Resolution for sampling the PDF on a grid.
    """
    # 1) Build the curve as a Shapely LineString
    curve = LineString(coordinates)

    # 2) Create the buffer (shadow) around the curve
    shadow = curve.buffer(r)

    # 3) Compute the area
    area_of_region = shadow.area

    # 4) Define a uniform PDF over that region
    def pdf_function(x, y):
        # Return 1/area if (x,y) is inside the buffer, else 0
        pt = Point(x, y)
        return (1.0 / area_of_region) if shadow.contains(pt) else 0.0

    # --- Plot the region ---
    x_vals, y_vals = zip(*coordinates)
    shadow_x, shadow_y = shadow.exterior.xy

    fig, ax = plt.subplots(figsize=(10,8))
    ax.set_facecolor("#8F689A")  # dark purple background

    # Fill the polygon
    plt.fill(shadow_x, shadow_y, color='green', alpha=0.6,
              edgecolor='#4B0082', linewidth=2, label='Uniform region')

    # Plot the original curve in blue
    plt.plot(x_vals, y_vals, color='blue', linewidth=2, label='Curve (M,W)')

    plt.xlabel("M")
    plt.ylabel("W")
    plt.title(f"Uniform PDF around (M,W) curve, radius = {r}")
    plt.legend()
    plt.savefig('prior_uniform_2022.pdf', format='pdf', bbox_inches='tight', dpi=1000)
    plt.show()

    # --- Sample the PDF on a grid ---
    minx, miny, maxx, maxy = shadow.bounds
    xs = np.linspace(minx, maxx, grid_res)
    ys = np.linspace(miny, maxy, grid_res)
    X, Y = np.meshgrid(xs, ys)

    Z = np.zeros_like(X)
    for i in range(grid_res):
        for j in range(grid_res):
            Z[i, j] = pdf_function(X[i, j], Y[i, j])

    # --- Save to NPZ ---
    np.savez("2d_prior_test_2022_short.npz",
              area=area_of_region,
              X=X,
              Y=Y,
              Z=Z)
    #print("Saved '2d_prior_test_2022.npz_short' with keys: [area, X, Y, Z].")

    return area_of_region, shadow, pdf_function


#plot_curve_with_fixed_shadow_and_uniform_pdf(parameter, r=100, grid_res=100)
#area, shadow_polygon, pdf_func = plot_curve_with_fixed_shadow_and_uniform_pdf(parameter, r=100, grid_res=400)

file_path = 'infoupdated.csv'

T=[]
H=[]
P=[]
#for year in years:
year=2022
T=list(get_temperature_for_year(file_path, year)[1])
H=list(get_temperature_for_year(file_path, year)[2])
P=list(get_temperature_for_year(file_path, year)[3])

# plt.plot(T)
# plt.show()
cases=[]
























prediction_lead=7
from statsmodels.tsa.statespace.sarimax import SARIMAX


# def forecast_ar(data, steps=7):
#     """
#     Forecast using SARIMA model.

#     Parameters:
#         data: array-like
#             Time series data.
#         steps: int
#             Number of steps to forecast.
#         order: tuple
#             (p,d,q) parameters for non-seasonal ARIMA.
#         seasonal_order: tuple
#             (P,D,Q,s) parameters for seasonal component, s is seasonal periodicity.

#     Returns:
#         forecast: array-like
#             Forecasted values.
#     """
#     order=(1,1,1)
#     seasonal_order=(1,1,1,7)
#     model = SARIMAX(data, order=order, seasonal_order=seasonal_order,
#                     enforce_stationarity=False, enforce_invertibility=False)
#     results = model.fit(disp=False)
#     forecast = results.predict(start=len(data), end=len(data)+steps-1)
    
#     return forecast

def forecast_ar(data, steps=7, lags=7):
    model = AutoReg(data, lags=7).fit()
    forecast = model.predict(start=len(data), end=len(data) + steps - 1)
    return forecast

# for i in range(340,365):   
for i in range(starting_time,365,prediction_lead):
    #print("i",i-(3*365))
    # print("i",i)
    d=i
    inuse=T[0:i]
    # model = AutoReg(inuse, lags=prediction_lead)
    # model_fitted = model.fit()
    predictions = forecast_ar(inuse, steps=prediction_lead, lags=7)
    # plt.plot(predictions)
    # plt.show()
    TT=list(predictions)
    # dama=T[3*365:i]+TT
    dama=T[0:i]+TT
    C=meanC

    
    
    inuseH=H[0:i]
    # modelH = AutoReg(inuseH, lags=prediction_lead)
    # model_fittedH = modelH.fit()
    predictionsH = forecast_ar(inuseH, steps=prediction_lead, lags=7)
    HH=list(predictionsH)
    #Humidity=H[3*365:i]+HH
    Humidity=H[0:i]+HH
    
    inuseP=P[0:i]
    # modelP = AutoReg(inuseP, lags=prediction_lead)
    # model_fittedP = modelP.fit()
    predictionsP = forecast_ar(inuseP, steps=prediction_lead, lags=7)
    PP=list(predictionsP)
    #Precipitation=P[3*365:+i]+PP
    Precipitation=P[0:i]+PP
    
    
    Coematrix = np.zeros((8,i+prediction_lead )) 
    for time in range(0, i+prediction_lead):
        Coematrix[0, time] = gammaMP(time)
        Coematrix[1, time] = pLSP(time)
        Coematrix[2, time] = muMP(time)
        Coematrix[3, time] = BRP(time)
        Coematrix[4, time] = alphaEP(time)
        Coematrix[5, time] = bEP(time)
        Coematrix[6, time] = PDRP(time)
        Coematrix[7, time] = VCP(time)
    
    def BMPPT(x,time):
        time=int(divmod(time, i+1)[1])
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
    t=np.linspace(0,i+prediction_lead,i+prediction_lead-1)
    x0=(3000000,0,0,0,10,10,10,10,10,100,100,100,100,100,100)
    x=odeint(BMPPT,x0,t)
    M=x[:,6]
    mosquito_profiles=np.cumsum(M)

    Dama = np.array(dama)
    Humidity = np.array(Humidity)
    Precipitation = np.array(Precipitation)
    
    # Now you can perform arithmetic operations
    BR_dic = np.cumsum(0.5 * (Dama + Humidity + Precipitation))

    
    Curve2022 = []

    for j in range(0, len(mosquito_profiles)):
        #print(j)
        Curve2022.append((mosquito_profiles[j], BR_dic[j]))


    # Save to a text file in two columns format
    with open("Curve-2022.txt", "w") as file:
        for item in Curve2022:
            file.write(f"{item[0]}\t{item[1]}\n")
    #grid_x, grid_y, grid_z = kernel_frequency_surface(all_coordinates, bw=accuracy)
    parameter=np.loadtxt('Curve-2022.txt')
    AAAA=parameter
    parameter=parameter[-1]
    #area, shadow_polygon, pdf_func = plot_curve_with_fixed_shadow_and_uniform_pdf(AAAA, r=250, grid_res=400)
    # print("AAAAAAAAAAAAAAAAAA",(list(parameter)[0],list(parameter)[1]))
    # print("K",(i-(3*365)))
    # if i % 7==0 and int(i/7)<52:
    #     k=i
    #     print(k/7)
    #     if Time[int(k/7)]>0 and int(k/7) <52 :
    #         print("updating time",k/7,Time[int(k/7)],k)
    #         # print(Time[int(k/7)],k)
    #         v=[(list(parameter)[0],list(parameter)[1])]#*(int(Time[int(k/7)]))
    #         # print("v",v)
            
    #         all_coordinates=all_coordinates+v 
    #         grid_x, grid_y, grid_z = kernel_frequency_surface(all_coordinates, bw=accuracy)


    #=========================================post==========================================
    #=========================================post==========================================
    #=========================================post==========================================
    #=========================================post==========================================
    #=========================================post==========================================
    
    def get_from_pdf(pdf_file, m, w):
        """
        Interpolate density values from the saved 2D PDF file, ensuring m and w are within bounds.
        """
        data = np.load(pdf_file)
        X = data['X']
        Y = data['Y']
        Z = data['Z']
    
        x_unique = np.unique(X)
        y_unique = np.unique(Y)
    
        if Z.shape == (len(y_unique), len(x_unique)):
            Z = Z.T
    
        # Clamp m and w to valid ranges
        m = np.clip(m, x_unique.min(), x_unique.max())
        w = np.clip(w, y_unique.min(), y_unique.max())
    
        interpolator = RegularGridInterpolator((x_unique, y_unique), Z)
        point = np.array([[m, w]])
        return interpolator(point)[0],x_unique.max(),y_unique.max(),x_unique.min(),y_unique.min()
    
    # pdf_file = "2d_pdf.npz"
    # print(get_from_pdf(pdf_file, m=5e12, w=4000))
    m=5e12
    w=4000
    
    #------------------------------------------------------------------------------
    #--------defining the posterior we and plotting posterior for given X=x--------
    #------------------------------------------------------------------------------
    
    # Posterior function
    
   
    
    year = 2022
    
    # ------------------------------------------------------------------
    # 1) Vectorized version of get_from_pdf:
    #    This function should load a 2D grid from `pdf_file` (X, Y, Z)
    #    and interpolate for all points in M_array, W_array simultaneously.
    #    You can use e.g. `scipy.interpolate.RegularGridInterpolator`
    #    or your own interpolation. Below is a skeleton you must adapt
    #    to match how `pdf_file` is actually stored.
    
    
    
    
    def get_from_pdf_array(pdf_file, M_array, W_array):
        """
        Vectorized version of get_from_pdf.
    
        Inputs:
          pdf_file (str): e.g. "2d_prior_test.npz"
          M_array, W_array: 2D arrays of shape (nrow, ncol)
                           specifying (m, w) points at which we want the PDF values.
    
        Returns:
          A 2D numpy array 'values' of shape (nrow, ncol),
          giving the PDF (or rate) at each (M_array[i,j], W_array[i,j]).
        """
        data = np.load(pdf_file)
        X = data['X']  # might be 1D or 2D
        Y = data['Y']
        Z = data['Z']  # the PDF or rate values
    
        # Convert X, Y to unique 1D arrays for interpolation.
        # (Same logic as in your 'get_from_pdf' function.)
        x_unique = np.unique(X)
        y_unique = np.unique(Y)
    
        # If Z has shape (len(y_unique), len(x_unique)), transpose it
        # so Z matches shape (len(x_unique), len(y_unique)).
        if Z.shape == (len(y_unique), len(x_unique)):
            Z = Z.T
    
        # Now Z should have shape (len(x_unique), len(y_unique))
        # and we build the interpolator with (x_unique, y_unique) in that order.
        interpolator = RegularGridInterpolator((x_unique, y_unique), Z)
    
        # Clamp the M_array and W_array to be within valid ranges
        M_clamped = np.clip(M_array, x_unique.min(), x_unique.max())
        W_clamped = np.clip(W_array, y_unique.min(), y_unique.max())
    
        # Build the list of points [(m, w), (m, w), ...] for all grid points
        points = np.column_stack((M_clamped.ravel(), W_clamped.ravel()))
    
        # Interpolate in one vectorized call
        vals_1d = interpolator(points)  # shape (#points,)
    
        # Reshape back to 2D
        values_2d = vals_1d.reshape(M_array.shape)
    
        return values_2d
    
    # ------------------------------------------------------------------
    # 2) Posterior function (vectorized)
    #    Instead of calling get_from_pdf(rate_file, m, w) for each single (m,w),
    #    we do it once for the entire 2D arrays.
    
    def posterior_array(M_array, W_array, x_value, prior_file="2d_prior_test_2022_short.npz", rate_file="2d_rate_2022_short.npz"):
        """
        Compute the posterior for entire 2D arrays M_array, W_array at once.
        
        posterior = PoissonPMF(x_value; mu=rate) * prior_pdf(m, w).
        """
        # 2.1) Get the 2D array of rate values
        r_values = get_from_pdf_array(rate_file, M_array, W_array)  # shape (nrow, ncol)
        
        # 2.2) Get the 2D array of prior pdf values
        prior_values = get_from_pdf_array(prior_file, M_array, W_array)  # shape (nrow, ncol)
        
        # 2.3) Poisson PMF for x_value, broadcast over entire r_values
        #      SciPy's poisson.pmf(k, mu=...) is vectorized for mu as an array
        pmf_values = poisson.pmf(x_value, mu=r_values)  # shape (nrow, ncol)
        
        # 2.4) posterior = pmf_values * prior_values
        return pmf_values #* prior_values


    # def posterior_array(M_array, W_array, x_value, prior_file="2d_prior_test_2022_short.npz", rate_file="2d_rate_2022_short.npz"):
    #     """
    #     Compute the posterior for entire 2D arrays M_array, W_array at once, ensuring:
    #       1. The posterior is zero outside the domain of the prior.
    #       2. The posterior is normalized to be a proper PDF.
        
    #     Posterior = PoissonPMF(x_value; mu=rate) * prior_pdf(m, w).
        
    #     Parameters:
    #         M_array, W_array: 2D arrays representing meshgrid of (m, w) values.
    #         x_value: Observed count for the Poisson likelihood.
    #         prior_file: Filename containing the prior PDF values.
    #         rate_file: Filename containing the Poisson rate values.
    
    #     Returns:
    #         2D array of the normalized posterior PDF.
    #     """
    #     # 1) Get the 2D array of rate values
    #     r_values = get_from_pdf_array(rate_file, M_array, W_array)  # shape (nrow, ncol)
        
    #     # 2) Get the 2D array of prior pdf values
    #     prior_values = get_from_pdf_array(prior_file, M_array, W_array)  # shape (nrow, ncol)
        
    #     # 3) Compute the Poisson PMF for x_value
    #     pmf_values = poisson.pmf(x_value, mu=r_values)  # shape (nrow, ncol)
        
    #     # 4) Compute raw posterior
    #     raw_posterior = pmf_values * prior_values  # shape (nrow, ncol)
        
    #     # 5) Ensure posterior has the same domain as the prior (zero outside)
    #     posterior_masked = np.where(prior_values > 0, raw_posterior, 0)
        
    #     # 6) Normalize the posterior to ensure it integrates to 1 (valid PDF)
    #     posterior_sum = np.sum(posterior_masked)
    #     if posterior_sum > 0:
    #         posterior_masked /= posterior_sum  # Normalize to make it a proper PDF
    
    #     return posterior_masked

    
    # ------------------------------------------------------------------
    # 3) Main code
    
    posterior_data = {}
    
    for i in range(tt, 32):
        pdf_file = "2d_prior_test_2022_short.npz"
    
        # We assume get_from_pdf() or get_from_pdf_array() can also return
        # the bounding box or min/max M,W, but let's just load it ourselves:
        data = np.load(pdf_file)
        # data should have X, Y, Z
        X_loaded = data["X"]  # Suppose shape is (ny, nx) or (nx,) if 1D
        Y_loaded = data["Y"]
        
        # If X_loaded and Y_loaded are 1D:
        maxm = X_loaded.max()
        minm = X_loaded.min()
        maxw = Y_loaded.max()
        minw = Y_loaded.min()
        
        #print("aaaaa", maxm)
    
        # Domain for M and W
        M = np.linspace(minm, maxm, 100)
        W = np.linspace(minw, maxw, 100)
        X, Y = np.meshgrid(M, W)  # shape (100,100)
    
        # 3.1) Vectorized posterior calculation
        x_value = i
        Z = posterior_array(X, Y, x_value,
                            prior_file="2d_prior_test_2022_short.npz",
                            rate_file="2d_rate_2022_short.npz")
    
        # 3.2) Normalize
        Z_sum = Z.sum()
        if Z_sum > 0:
            Z /= Z_sum
    
        # 3.3) Save posterior data as (m, w, density, x_value)
        m_values = X.ravel()
        w_values = Y.ravel()
        density_values = Z.ravel()
        x_values = np.full_like(m_values, fill_value=x_value, dtype=np.float64)
    
        posterior_array_4cols = np.column_stack((m_values, w_values, density_values, x_values))
        np.save(f'posterior_x_{x_value}_2022_short.npy', posterior_array_4cols)
        posterior_data[x_value] = posterior_array_4cols
    
        # 3.4) Plot
        # plt.figure(figsize=(10, 8))
        # contour = plt.contourf(X, Y, Z, levels=50, cmap='viridis', alpha=0.7)
        # cbar = plt.colorbar(contour)
        # cbar.set_label("Posterior Probability", fontsize=12)
    
        # Locate the maximum posterior
        max_pos = np.unravel_index(Z.argmax(), Z.shape)
        max_m = X[max_pos]
        max_w = Y[max_pos]
    
        #plt.plot(max_m, max_w, '*k', markersize=8, label=f'Max Posterior (x={x_value})')
        
        
        
        
        parameters = []
        parameters = np.loadtxt("Curve-2022.txt")
    
        param_x = [point[0] for point in parameters]
        param_y = [point[1] for point in parameters]
        # plt.scatter(param_x, param_y, color='black', alpha=0.5, s=12, label='Parameters')
        # plt.text(
        #     max_m, max_w, f'{x_value}',
        #     color='white', fontsize=10,
        #     ha='left', va='bottom',
        #     bbox=dict(boxstyle="round,pad=0.3", fc="black", ec="none", alpha=0.7)
        # )
    
        # plt.title("Posterior Probability Surface with Parameters", fontsize=14)
        # plt.xlabel("Mosquito Axis (M)", fontsize=12)
        # plt.ylabel("Weather Axis (W)", fontsize=12)
        # plt.legend(fontsize=10)
        # plt.tight_layout()
        # plt.show()
    



    
    #------------------------------------------------------------------------------
    #------------------------------------------------------------------------------
    #--------------------------chwck poaterior fitted
    
    def preload_posteriors(x_values=range(tt,32)):
        """
        Preloads all posterior_x_{x_val}.npy files, reshapes them into
        a 2D array, and builds a RegularGridInterpolator for each x_val.
        
        Returns a dictionary:
          interpolators[x_val] = (f, min_m, max_m, min_w, max_w)
        where:
          f is a RegularGridInterpolator object
          min_m, max_m, min_w, max_w are the bounding box to clamp queries.
        """
        interpolators = {}
    
        for x_val in x_values:
            filename = f"posterior_x_{x_val}_2022_short.npy"
            if not os.path.exists(filename):
                raise FileNotFoundError(f"File '{filename}' not found.")
    
            data = np.load(filename)
            # data shape: (N,4): columns = [m_i, w_i, density_i, x_i]
    
            all_m = data[:, 0]
            all_w = data[:, 1]
            all_d = data[:, 2]
    
            # Identify the unique sorted M and W
            unique_m = np.unique(all_m)
            unique_w = np.unique(all_w)
            nM = len(unique_m)
            nW = len(unique_w)
    
            # Build a 2D array 'density_2d'
            density_2d = np.zeros((nW, nM), dtype=float)
    
            # Fill the array by matching each row in data to the correct [i_w, i_m]
            # This is typically O(N) = 10000 if it's a 100x100 grid, but done once.
            for row in data:
                m_i, w_i, dens_i, _ = row
                i_m = np.searchsorted(unique_m, m_i)
                i_w = np.searchsorted(unique_w, w_i)
                density_2d[i_w, i_m] = dens_i
    
            # Build the RegularGridInterpolator
            # Note the domain order is (w, m) if we interpret density_2d as [i_w, i_m].
            f = RegularGridInterpolator(
                (unique_w, unique_m),    # (y-axis array, x-axis array)
                density_2d,              # shape (nW, nM)
                bounds_error=False,
                fill_value=0.0           # Extrapolate => 0
            )
    
            # Store in dictionary
            min_m, max_m = unique_m.min(), unique_m.max()
            min_w, max_w = unique_w.min(), unique_w.max()
            interpolators[x_val] = (f, min_m, max_m, min_w, max_w)
    
        return interpolators
    
    
    
    
    def find_best_x_for_point(m, w, interpolators):
        """
        Given (m, w) and a dictionary of interpolators built by 'preload_posteriors',
        returns the x_val with the highest interpolated density at (m, w).
        
        Args:
          m, w (floats): the query point
          interpolators (dict): as returned by preload_posteriors
        
        Returns:
          best_x  (int): The x-value with the highest density
          best_density (float): The density at that (m, w)
        """
        best_x = None
        best_density = -1.0
    
        # Loop over all x_vals we cached
        for x_val, (f, min_m, max_m, min_w, max_w) in interpolators.items():
            # Clamp (m, w) to valid range
            mm = np.clip(m, min_m, max_m)
            ww = np.clip(w, min_w, max_w)
    
            # Evaluate the interpolator
            density_at_mw = f((ww, mm))  # shape=() single float
    
            if density_at_mw > best_density:
                best_density = density_at_mw
                best_x = x_val
    
        return best_x, float(best_density)
    
    
    
    
    
    
    interpolators_dict = preload_posteriors(x_values=range(tt,32))
    
    
    
    
    
    
    
    
    
    
    
    #file_path = 'Weather_Mosquito_Profiles_2022.txt'
    #file_path = 'Curve2022_short.txt'
    # loaded_data = parameter#np.loadtxt(file_path, skiprows=1) [-1] # Skip the header row
    
    # Extract mosquito and weather profiles
    # loaded_mosquito_profiles = loaded_data[:, 0]
    # loaded_weather_profiles = loaded_data[:, 1]
    
    # print(f"Loaded Mosquito Profiles for {year}: {loaded_mosquito_profiles}")
    # print(f"Loaded Weather Profiles for {year}: {loaded_weather_profiles}")
    # Example usage
    x_values = list(range(tt, 32))  # The range of x_values for which posteriors are saved
    # m0 = 5.5e12  # Example m-coordinate
    # w0 = 4000    # Example w-coordinate
    
    # max_x_value, max_density = find_highest_posterior(x_values, m0, w0)
    # print(f"The highest posterior density at ({m0}, {w0}) is for x_value={max_x_value}, with density={max_density}.")

    #for i in range(starting,starting_time+1):
    m=parameter[0]#loaded_mosquito_profiles[i]
    # print("m",m)
    w=parameter[1]#loaded_weather_profiles[i]
    #print(m,w)
    # print(find_best_x_for_point(m, w, interpolators_dict)[0])
    # cases.append(find_best_x_for_point(m, w, interpolators_dict)[0])
    print("day",d+prediction_lead,"prediction",find_best_x_for_point(m, w, interpolators_dict)[0])
    cases.append(find_best_x_for_point(m, w, interpolators_dict)[0])
    if (d) % 7==0 and int((d)/7)<52:
        k=d
        #print("updated to ",k/7)
        if Time[int(k/7)]>0 and int(k/7) <52 :
            r=int(k)
            print("updating day",k,"week",int(k/7),"weight",Time[int(k/7)])
            # print(Time[int(k/7)],k)
            v=[(list(parameter)[0],list(parameter)[1])]#*(int(Time[int(k/7)]))
            # print("v",v)
            
            all_coordinates=all_coordinates+v 
            grid_x, grid_y, grid_z = kernel_frequency_surface(all_coordinates, bw=accuracy)

    
     
    
    
    
print(cases)    
 



















def visualize_risk(predicted_cases, mosquito_profile,year, A, B):
    """
    Visualize the risk of West Nile disease infection in Orange County with mosquito profiles.
    
    This function creates a figure with two subplots:
    
      (1) The top subplot shows:
          - Predicted daily cases as color-coded bars (with colors ranging from yellow to red)
          - Weekly data (Time) as black stars (assumed one per week, plotted at x = 7 * week_index)
          - The mosquito profile as a blue line on a secondary y-axis.
          
      (2) The bottom subplot shows, for each provided Poisson parameter a (from A):
          - A vertical green line starting at 0 and ending at the 99th percentile of Poisson(a)
          - A red circle marking the corresponding b value (from B) on that line.
          - Each vertical line is positioned at x = 7 * week_index to align with the weekly observations.
          
    Parameters:
        predicted_cases (list or array): Predicted daily case counts for days [starting, ending).
        mosquito_profile (list or array): Mosquito profile values (one per day).
        Time (list or array): Weekly observed values (assumed to have 52 elements).
        year (int): The year for which the visualization is generated.
        A (list or array): A set of a–values (used as the parameter for a Poisson PMF).
        B (list or array): A set of b–values to be marked on the corresponding vertical lines.
        starting (int): Starting day index (default = 0).
        ending (int): Ending day index (exclusive; default = 364).
        
    Returns:
        None
    """

    # Prepare the x-axis for days in the top plot
    days = np.arange(starting, 365)  # e.g. 0, 1, 2, ..., 363


    # Create a colormap and normalization for predicted cases (assumes cases from 1 to 32)
    cmap = plt.colormaps.get_cmap('YlOrRd')  # Yellow-to-red colormap
    norm = Normalize(vmin=1, vmax=32)

    # Create a figure with two subplots (using gridspec so the bottom plot is smaller)
    fig = plt.figure(figsize=(13, 12))
    gs = fig.add_gridspec(2, 1, height_ratios=[3, 1])
    ax1 = fig.add_subplot(gs[0, 0])  # Top subplot: main plot
    #ax3 = fig.add_subplot(gs[1, 0])  # Bottom subplot: Poisson vertical lines

    # ----- TOP PLOT: Original Visualization -----
    # 1) Plot predicted daily cases as color-coded bars.
    max_Time = max(Time) if len(Time) > 0 else 1
    for ii in range(0,starting-1):
        ax1.bar(ii,height=max_Time, width=1, color='green')
        
        
    for day, cases in zip(days, predicted_cases[:365]):
        bar_color = cmap(norm(cases))
        ax1.bar(day, height=max_Time, width=1, color=bar_color,
                align='edge', edgecolor='none')
    
    # 2) Add a horizontal colorbar for the predicted cases.
    sm = ScalarMappable(cmap=cmap, norm=norm)
    cbar = fig.colorbar(sm, ax=ax1, orientation='horizontal', fraction=0.05, pad=0.2)
    cbar.set_label("Severity (Cases)", fontsize=12)
    
    # 3) Plot the weekly observed cases as black stars.
    weeks = np.arange(len(Time))
    star_x = weeks * 7  # Each week starts at day = 7 * week_index
    ax1.plot(star_x, Time, '*k', label='Weekly Observed Cases')
    
    # 4) Add the mosquito profile on a second y-axis.
    ax2 = ax1.twinx()
    ax2.plot(np.arange(0, 364), mosquito_profile, color='blue', linewidth=1.5,
              label="Mosquito Profile")
    ax2.set_ylabel("Mosquito Profile", fontsize=12, color='blue')
    ax2.tick_params(axis='y', labelcolor='blue')
    
    # 5) Set title, x–axis label, and ticks for the top plot.
    ax1.set_title(f"Short-term Risk Prediction for Severity of WND, Orange County with Mosquito Profile ({year})",
                  fontsize=16)
    ax1.set_xlabel("Day of the Year", fontsize=12)
    ax1.set_xticks(np.arange(starting, 365, step=30))
    
    # 6) Combine legends from ax1 and ax2.
    lines_labels_1 = ax1.get_legend_handles_labels()
    lines_labels_2 = ax2.get_legend_handles_labels()
    lines = lines_labels_1[0] + lines_labels_2[0]
    labels = lines_labels_1[1] + lines_labels_2[1]
    ax1.legend(lines, labels, loc='upper left')
    
  
    plt.tight_layout()
    plt.savefig(f'severity-short{year}.pdf',dpi=1000)
    plt.show()



# def calculate_log_likelihood(A, B):
#     """
#     Compute the log-likelihood score by evaluating how well the observations (B) fit into
#     the expected Poisson distribution defined by parameters (A).
    
#     Fixes:
#     - If a_val == 0 and b_val == 0, log-likelihood = 0 (since P(X=0) = 1).
#     - If b_val is outside the 99% domain, assign a log-likelihood of -10.
#     - Prevents log(0) by ensuring probabilities are valid.

#     Parameters:
#         A (list or array): Poisson parameters for each of the 52 weeks.
#         B (list or array): Observed values corresponding to each week.

#     Returns:
#         total_log_likelihood (float): The sum of all log(P) values.
#     """
#     log_likelihoods = []

#     for week_index, (a_val, b_val) in enumerate(zip(A, B)):
#         # Compute 99% quantile for Poisson(a_val)
#         quantile_99 = int(poisson.ppf(0.99, a_val))  # Ensure integer for domain limit
        
#         if b_val >= 0:
#             # Handle special case: a_val == 0 and b_val == 0
#             if a_val == 0 and b_val == 0:
#                 log_likelihoods.append(0)  # log(1) = 0
#                 print(f"Week {week_index + 1}: Log-likelihood = 0 (a_val = 0, b_val = 0)")
#                 continue
#             if a_val != 0 and b_val == 0:
#                 log_likelihoods.append(-10)  
#                 print(f"Week {week_index + 1}: Log-likelihood = -10 (a_val = {a_val}, b_val = 0)")
#                 continue
    
#             # If b_val is outside the domain (greater than 99% quantile), assign log-likelihood = -10
#             if b_val > quantile_99:
#                 log_likelihoods.append(-10)
#                 print(f"Week {week_index + 1}: Log-likelihood = -10 (b_val > 99% quantile)")
#                 continue
    
#             # Define bin edges using Poisson CDF
#             bin_edges = [int(poisson.ppf(p, a_val)) for p in np.linspace(0.01, 0.99, 10)]
            
#             # Ensure unique and sorted bin edges
#             bin_edges = sorted(set(bin_edges))
    
#             # Find the bin in which b_val falls
#             bin_index = None
#             for i in range(len(bin_edges) - 1):
#                 if bin_edges[i] <= b_val < bin_edges[i + 1]:
#                     bin_index = i
#                     break
    
#             # If b_val is smaller than the first bin, assign it to the first bin
#             if bin_index is None:
#                 bin_index = 0 if b_val < bin_edges[0] else len(bin_edges) - 2
    
#             # Compute probability for the bin
#             P_bin = poisson.cdf(bin_edges[bin_index + 1], a_val) - poisson.cdf(bin_edges[bin_index], a_val)
    
#             # Avoid log(0) errors but check for perfect probability
#             log_likelihood = np.log(max(P_bin, 1e-10)) if not np.isclose(P_bin, 1.0) else 0
            
#             log_likelihoods.append(log_likelihood)
#             print(f"Week {week_index + 1}: Log-likelihood = {log_likelihood:.6f}")

#     # Sum all log(P) values
#     total_log_likelihood = sum(log_likelihoods)/52
    
#     print(f"\nTotal Log-Likelihood: {total_log_likelihood:.6f}")
    
#     return total_log_likelihood


def calculate_log_likelihood(A, B):
    """
    Compute the log-likelihood score by evaluating how well the observations (B) fit into
    the expected Poisson distribution defined by parameters (A), now aggregated monthly.

    Fixes:
    - Converts weekly data (52 weeks) into **monthly totals**.
    - If a_val == 0 and b_val == 0, log-likelihood = 0 (since P(X=0) = 1).
    - If b_val is outside the 99% Poisson domain, assigns log-likelihood = -10.
    - Prevents log(0) errors by ensuring probabilities are valid.

    Parameters:
        A (list or array): Poisson parameters for each of the 52 weeks.
        B (list or array): Observed values corresponding to each week.

    Returns:
        total_log_likelihood (float): The mean log-likelihood per month.
    """

    # === Step 1: Convert Weekly Data to Monthly Data ===
    if len(A) != 52 or len(B) != 52:
        raise ValueError("Expected 52 weeks of data for A and B.")

    # Define the number of weeks per month (approximated)
    month_week_mapping = {
        1: 4, 2: 4, 3: 5, 4: 4, 5: 4, 6: 5,
        7: 4, 8: 4, 9: 5, 10: 4, 11: 4, 12: 5
        # 1: 4, 2: 4, 3: 4, 4: 4, 5: 4, 6: 4,
        # 7: 4, 8: 4, 9: 4, 10: 4, 11: 4, 12: 4,13:4
    }

    A_monthly = []
    B_monthly = []

    week_index = 0
    for month, num_weeks in month_week_mapping.items():
        A_monthly.append(sum(A[week_index: week_index + num_weeks]))  # Sum Poisson rates
        B_monthly.append(sum(B[week_index: week_index + num_weeks]))  # Sum observed cases
        week_index += num_weeks

    # === Step 2: Compute Log-Likelihood Score for Monthly Data ===
    log_likelihoods = []

    for month_index, (a_val, b_val) in enumerate(zip(A_monthly, B_monthly), start=1):
        # Compute the 99% quantile for Poisson(a_val)
        quantile_99 = int(poisson.ppf(0.99, a_val))  # Ensure integer for domain limit
        
        if b_val >= 0:
            # Handle special case: a_val == 0 and b_val == 0
            if a_val == 0 and b_val == 0:
                log_likelihoods.append(0)  # log(1) = 0
                print(f"Month {month_index}: Log-likelihood = 0 (a_val = 0, b_val = 0)")
                continue
            if a_val != 0 and b_val == 0:
                log_likelihoods.append(-10)  
                print(f"Month {month_index}: Log-likelihood = -10 (a_val = {a_val}, b_val = 0)")
                continue
    
            # If b_val is outside the domain (greater than 99% quantile), assign log-likelihood = -10
            if b_val > quantile_99:
                log_likelihoods.append(-10)
                print(f"Month {month_index}: Log-likelihood = -10 (b_val > 99% quantile)")
                continue
    
            # Define bin edges using Poisson CDF
            bin_edges = [int(poisson.ppf(p, a_val)) for p in np.linspace(0.01, 0.99, 10)]
            
            # Ensure unique and sorted bin edges
            bin_edges = sorted(set(bin_edges))
    
            # Find the bin in which b_val falls
            bin_index = None
            for i in range(len(bin_edges) - 1):
                if bin_edges[i] <= b_val < bin_edges[i + 1]:
                    bin_index = i
                    break
    
            # If b_val is smaller than the first bin, assign it to the first bin
            if bin_index is None:
                bin_index = 0 if b_val < bin_edges[0] else len(bin_edges) - 2
    
            # Compute probability for the bin
            P_bin = poisson.cdf(bin_edges[bin_index + 1], a_val) - poisson.cdf(bin_edges[bin_index], a_val)
    
            # Avoid log(0) errors but check for perfect probability
            log_likelihood = np.log(max(P_bin, 1e-10)) if not np.isclose(P_bin, 1.0) else 0
            
            log_likelihoods.append(log_likelihood)
            print(f"Month {month_index}: Log-likelihood = {log_likelihood:.6f}")

    # Compute the mean log-likelihood across all months

    total_log_likelihood = sum(log_likelihoods) / len(log_likelihoods)
    
    print(f"\nTotal Monthly Log-Likelihood: {total_log_likelihood:.6f}")
    
    return total_log_likelihood

   
MM=list(np.loadtxt("./MMosquito_Profile_2022.txt"))

#visualize_risk(cases, MM, 2022)

def weekly_mean(data):
    return [np.mean(data[i:i+7]) for i in range(0, len(data), 7)]



def moving_average(X, window_size=3):
    if window_size < 1:
        raise ValueError("Window size must be at least 1")
    
    # Extend the array at the beginning and end to maintain length
    pad_size = window_size // 2
    X_padded = np.pad(X, (pad_size, pad_size), mode='edge')
    
    # Compute the moving average using convolution
    kernel = np.ones(window_size) / window_size
    smoothed = np.convolve(X_padded, kernel, mode='valid')
    
    return smoothed.tolist()

B=Time



A2=cases

A1=[0 for i in range(0,int(starting/7))]
#A1=[0 for i in range(0,31)]
#A2=[cases[i] for i in range(0,len(cases)) if i % 7 == 0]  
#A2=[1, 3, 2, 7, 9, 3, 4, 7, 5, 8, 13, 18, 8, 14, 16, 3, 2, 1, 0, 0, 0]

A=[0] * (52-len(A2)) + A2#[1, 3, 2, 7, 9, 3, 4, 7, 5, 8, 13, 18, 8, 14, 16, 3, 2, 1, 0, 0]

print(len(A))



AA=[0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 1, 1, 2, 2, 2, 3, 3, 3, 3, 2, 2, 2, 1, 1, 1, 0, 0, 0, 0]
# print(calculate_log_likelihood(AA[int(starting/7):], B[int(starting/7):]))
# y=calculate_log_likelihood(AA[int(starting/7):], B[int(starting/7):])

#print(len(A),len(AA))


x=calculate_log_likelihood(A, B)
y=calculate_log_likelihood(AA, B)



def plot_model_comparison(A, AA, B, year, x, y):
    """
    Plots the given model results and reported values for a given year.
    
    - A: Results of the probabilistic model (line plot).
    - AA: Results of the negative binomial model (line plot).
    - B: Reported real values (bar plot, black, thicker).
    - year: The year for labeling the plot.
    - x: Score of the probabilistic model.
    - y: Score of the negative binomial model.
    
    Saves the plot as a PDF file named 'model_comparison_YEAR.pdf'.
    """
    if len(A) != len(AA) or len(A) != len(B):
        raise ValueError("All input lists must have the same length.")

    weeks = np.arange(1, len(A) + 1)

    # Create the plot
    plt.figure(figsize=(8, 5))
    
    # Plot real cases as black bars (thicker)
    plt.bar(weeks, B, color='black', alpha=0.7, label='Reported Values', width=0.6)

    # Plot model results as lines with scores in the legend
    plt.plot(weeks, A, marker='o', linestyle='-', label=f'Probabilistic Model Result (Score = {np.round(x,decimals=2)})', linewidth=2)
    plt.plot(weeks, AA, marker='s', linestyle='--', label=f'Negative Binomial Model (Score = {np.round(y,decimals=2)})', linewidth=2)
    # Labels and title
    plt.xlabel('Week',fontsize=14)
    plt.ylabel('Number of Cases',fontsize=14)
    plt.title(f'Orange County-Californoia-{year}', fontsize=15)
    plt.legend()
    
    # Save as PDF
    filename = f"Short-term_model_comparison_oneweek{year}.pdf"
    plt.savefig(filename, format="pdf")

    # Show the plot
    plt.show()

    print(f"Plot saved as {filename}")
year=2022
plot_model_comparison(A, AA, B, year,x,y)