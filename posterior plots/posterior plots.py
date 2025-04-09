# -*- coding: utf-8 -*-
"""
Created on Sun Jan 26 19:18:01 2025

@author: shosseini
"""

import numpy as np
import matplotlib.pyplot as plt
from scipy.integrate import odeint
#from matplotlib import pyplot as plt
from scipy.stats import gamma
from scipy.stats import gaussian_kde
import pandas as pd
from scipy.optimize import minimize
from scipy.stats import gamma, norm, multivariate_normal, chi2, kstest
from scipy.interpolate import RegularGridInterpolator

from matplotlib.colors import Normalize
from matplotlib.cm import ScalarMappable
import matplotlib.colors as mcolors
from scipy.stats import poisson
from scipy.interpolate import RectBivariateSpline
from scipy.interpolate import griddata

#---------------------------------------------------------------------------------
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
    
    return result


# File path and years to process
file_path = 'info.csv'





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


Coematrix=np.zeros([9,365])   


weeks_total_cases = np.zeros(52)  # Assuming 52 weeks in a year
years = [year for year in range(2006, 2021) if year != 2010]
specific_year = 2014  # Year to compare
specific_year_cases = None  # Placeholder for the specific year's data
for year in years:
    try:
        case_data = np.loadtxt(f"./orangecases{year}.txt")
        if year == specific_year:
            specific_year_cases = case_data  # Save the specific year's data
        for week, cases in enumerate(case_data):
            weeks_total_cases[week] += cases  # Accumulate cases per week
    except FileNotFoundError:
        print(f"Data for {year} not found. Skipping this year.")

years = list(range(2006, 2021))  # Years from 2006 to 2020 (excluding 2021)

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


# for num_cases, occurrences in sorted(cases_occurrences.items()):
#     print(f"Cases: {num_cases}, Occurrences: {occurrences}")

def moving_average(data, window_size):
    """
    Calculate the moving average of a list.

    Parameters:
        data (list or array-like): Input list of numbers.
        window_size (int): Size of the moving window.

    Returns:
        list: Moving average values.
    """
    if not data or window_size <= 0:
        raise ValueError("Data must be non-empty and window_size must be a positive integer.")
    
    if window_size > len(data):
        raise ValueError("Window size must not exceed the length of the data.")
    
    # Calculate the moving average
    moving_avg = []
    for i in range(len(data) - window_size + 1):
        window = data[i:i + window_size]
        moving_avg.append(sum(window) / window_size)
    
    return moving_avg

#--------------------------------------------------------------------------------


years = [year for year in range(2006, 2021) if year != 2010]
# Dictionary to store cumulative mosquito profiles (M) by year
mosquito_profiles = {}
Weather={}
mos={}


        
 
#years=[2007]
    

for year in years:
    print(year)
    dama = list(np.loadtxt(f"./T{year}.txt"))
    H = list(np.loadtxt(f"./H{year}.txt"))
    C = list(np.loadtxt(f"./Cf{year}.txt"))
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

    # Integrate BMPPT to compute mosquito profiles
    t = np.linspace(0, 365, 364)
    x0 = (3000000, 0, 0, 0, 10, 10, 10, 10, 10, 100, 100, 100, 100, 100, 100)
    x = odeint(BMPPT, x0, t)
    incB=1/3
    rC=1/6
    muwnvB=0.9
    muB=(0.15/365)

    # Compute cumulative mosquito values (M)
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
    
    mosquito_profiles[year] =M
    Weather[year]=np.cumsum(get_temperature_for_year(file_path, year))
    mos[year]=Mos
    # temperatures = get_temperature_for_year(file_path, year)
    curve=[]
    for i in range(0,364):
        curve.append((mosquito_profiles[year][i], Weather[year][i]))
  
#--------------------------------------------------------------------    
#--------------------------------------------------------------------------------
    # Define the years to process
    # Define years to process
    # years = [year for year in range(2006, 2021) if year != 2010]
    # print("Number of years processed:", len(years))
    
    # Initialize Cases array (14 rows for 2006–2020 excluding 2010, 365 days)
    Cases = np.zeros((14, 365))
    
    # Populate Cases with data
    for y in years:
        Casetime = list(np.loadtxt(f"./orangecases{y}.txt"))
        for i in range(0, 52):
            Cases[years.index(y), i * 7] = Casetime[i]  # Use the index of `year` in the `years` list
    
    # Define function to get case info for a specific year
    
    def case_number_byyear(time_year):
        """
        Retrieve daily cases for a specific year.
    
        :param time_year: Year to retrieve data for (e.g., 2014)
        :return: Daily cases for the specified year
        """
        if time_year not in years:
            raise ValueError(f"Year {time_year} is not in the list of processed years: {years}")
        # Use the index of the year in the `years` list
        year_index = years.index(time_year)
        return Cases[year_index, :]
    
    
    
    
    
    #------------------------------------------------------------------------------
    #----------------------------------------reading pdf from  another file
    #------------------------------------------------------------------------------
    #Load the saved PDF
    #------------------------------------------2-d prior---------------------------
    #------------------------------------------------------------------------------
    #-----------------------------------shows tha countor of 99 percent---------
    
    
    
    
    
    
    
    # Load the 2D PDF data
    data = np.load("2d_pdf.npz")  # File of the 2D PDF
    X = data['X']
    Y = data['Y']
    Z = data['Z']
    
    # Ensure X and Y are flattened for min/max operations
    X_flat = X.flatten()
    Y_flat = Y.flatten()
    
    # Plot the PDF
    # plt.figure(figsize=(10, 8))
    # contour = plt.contourf(X, Y, Z, levels=50, cmap='viridis', alpha=0.7)
    # plt.colorbar(contour, label='Fitted PDF Density')
    
    # Calculate the total density
    # total_density = np.sum(Z)
    
    # # Find the contour level enclosing 99% of the density
    # sorted_Z = np.sort(Z.flatten())[::-1]  # Flatten and sort in descending order
    # cumulative_density = np.cumsum(sorted_Z) / total_density
    # level_99 = sorted_Z[np.argmax(cumulative_density >= 0.99)]  # Density value at 99%
    
    # # Plot the 99% contour
    # contour_99 = plt.contour(X, Y, Z, levels=[level_99], colors='red', linewidths=2, linestyles='--')
    # contour_99.collections[0].set_label('99% Contour')  # Manually set the label for legend
    
    # # Read "Curve2021.txt" and plot it
    # curve_data =np.array(curve)# np.loadtxt(f"Curve{year}.txt")  # Load curve data (assuming two columns: X and Y)
    # curve_x, curve_y = curve_data[:, 0], curve_data[:, 1]
    
    
    
    
    
    # plt.plot(curve_x, curve_y, color='blue', linewidth=2, label=f'Curve from Curve{year}.txt')
    
    # # Add labels, title, and legend
    # plt.xlabel('X')
    # plt.ylabel('Y')
    # plt.title(f'Contour Enclosing 99% of the PDF with data for {year}')
    # plt.legend(fontsize=10)
    # plt.tight_layout()
    
    # # Use flattened X and Y for setting limits
    # plt.xlim(X_flat.min(), X_flat.max())
    # plt.ylim(Y_flat.min(), Y_flat.max())
    
    # plt.show()
    
    
    #------------------------------------------------------------------------------
    #-------function to check if the chosen points are insiude the countor (domain)
    #-------------------------------of 99 percent or not---------------------------

    
    # def is_inside_contour(m, w):
    #     # Clamp (m, w) to valid ranges
    #     m = np.clip(m, X.min(), X.max())
    #     w = np.clip(w, Y.min(), Y.max())
        
    #     # Interpolate Z at (m, w)
    #     interp_func = RectBivariateSpline(X[0, :], Y[:, 0], Z)
    #     z_value = interp_func(m, w)[0, 0]  # Interpolated Z value
        
    #     # Check if the Z value is greater than or equal to the 99% level
    #     return 1 if z_value >= level_99 else 0
    
    
    
    def get_from_pdf(pdf_file, m, w):
        """
        Load a precomputed 2D rate from an .npz file and interpolate
        to find the rate at point (m, w).
    
        Parameters:
            pdf_file (str): Path to the .npz file containing arrays X, Y, Z.
            m (float): Mosquito value.
            w (float): Weather value.
    
        Returns:
            float: The rate at (m, w).
        """
        # Load 2D rate data
        data = np.load(pdf_file)
        X = data['X']
        Y = data['Y']
        Z = data['Z']
    
        # Ensure unique sorted 1D arrays for x_unique and y_unique
        x_unique = np.unique(X)
        y_unique = np.unique(Y)
    
        # If Z is shaped (len(y_unique), len(x_unique)), transpose it
        if Z.shape == (len(y_unique), len(x_unique)):
            Z = Z.T
    
        # Create the interpolator
        interpolator = RegularGridInterpolator((x_unique, y_unique), Z)
    
        # Evaluate the rate at the single point (m, w)
        point = np.array([[m, w]])
        rate = interpolator(point)[0]
        return rate
    
    
    pdf_file = "2d_pdf.npz"
    
    # # Example (m, w) values
    m = 5.209648826723815e12
    w = 3321.8949999999995
    # print("pdf",get_from_pdf(pdf_file, m, w))
    #------------------------------------------------------------------------------
    # -------------------plots the file and show a point on that read file---------
    #-----------------------------------------------------------------------------
    
    def plot_pdf_with_point(pdf_file, m, w):
        """
        Load the rate, compute the rate at (m, w), and plot the rate contour
        with the point (m, w) marked, along with its value.
    
        Parameters:
            pdf_file (str): Path to the .npz file containing arrays X, Y, Z.
            m (float): Mosquito value.
            w (float): Weather value.
        """
        # 1) Load data
        data = np.load(pdf_file)
        X = data['X']
        Y = data['Y']
        Z = data['Z']
    
        # 2) Create unique sorted axes
        x_unique = np.unique(X)
        y_unique = np.unique(Y)
    
        # If Z is shaped (len(y_unique), len(x_unique)), transpose it
        if Z.shape == (len(y_unique), len(x_unique)):
            Z = Z.T
    
        # 3) Prepare meshgrid for plotting
        Xg, Yg = np.meshgrid(x_unique, y_unique, indexing='ij')
    
        # 4) Plot the rate as a contour
        plt.figure(figsize=(10, 8))
        contour = plt.contourf(Xg, Yg, Z, levels=50, cmap='viridis')
        cb = plt.colorbar(contour)
        cb.set_label('Rate')
    
        # 5) Compute the rate at (m, w)
        rate = get_from_pdf(pdf_file, m, w)
    
        # 6) Mark the point on the plot
        #plt.plot(m, w, 'ro', markersize=8, label=f'(m, w) Point: Rate = {int(rate)}')
    
    
        # 7) Annotate the rate value near the point
        plt.text(
            m, w, f'{int(rate)}',
            color='white', fontsize=10,
            ha='left', va='bottom',
            bbox=dict(boxstyle="round,pad=0.3", fc="black", ec="none", alpha=0.7)
        )
    
        # 8) Final plot labeling
        plt.xlabel("Mosquito Axis (m)")
        plt.ylabel("Weather Axis (w)")
        plt.title("2D Rate Contour with (m, w) Marked")
        plt.legend()
        plt.xlim(1.9e12,7e12)
        plt.ylim(2000,5000)
        plt.tight_layout()
        plt.show()
    
    pdf_file = "2d_rate.npz"
    #pdf_file = "2d_pdf.npz"
    m_value = 4.209648826723815e12
    w_value = 3321.8949999999995
    plot_pdf_with_point(pdf_file, m_value, w_value)
    
    
    
    
    #------------------------------------------------------------------------------
    #------------------------------------------------------------------------------
    #                    Function to get density from the PDF file
    # The `get_from_pdf` function interpolates the density value of a 2D probability
    #                     distribution at a given point (m, w) 
    # using a pre-saved PDF grid. It ensures the input coordinates are within bounds
    #                     and returns the interpolated density 
    #                     along with the grid's maximum m and w values.
    #------------------------------------------------------------------------------
    
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
        return interpolator(point)[0],x_unique.max(),y_unique.max()
    
    pdf_file = "2d_pdf.npz"
    print(get_from_pdf(pdf_file, m=5e12, w=4000))
    

#------------------------------------------------------------------------------
#--------defining the posterior we and plotting posterior for given X=x--------
#------------------------------------------------------------------------------








    def is_point_inside_strip(m, w, year, tolerance=0.1):
        """
        Check if a given (m, w) point falls inside the uniform strip.
    
        Parameters:
            m (float): The m-coordinate of the point.
            w (float): The w-coordinate of the point.
            year (int or str): Year or identifier used in the saved strip file.
            tolerance (float): Maximum allowed distance to consider the point inside.
    
        Returns:
            bool: True if the point is inside the strip, False otherwise.
        """
        # Load the strip data
        data = np.load(f'2d_prior_pdf-{year}.npz')
        strip_x = data['strip_x']
        strip_y = data['strip_y']
    
        # Compute distances from (m, w) to all points in the strip
        distances = np.sqrt((strip_x - m)**2 + (strip_y - w)**2)
    
        # Check if the minimum distance is within the tolerance
        return np.any(distances < tolerance)

    def posterior(m, w, x):
        """
        Compute the posterior value for given m, w, and x.
        """
        #pdf_file="2d_prior_pdf-mixed.npz"
        pdf_file=f"2d_prior_pdf-{year}.npz"
        #pdf_file = "2d_prior_pdf.npz"
        rate_file = "2d_rate.npz"
        r = get_from_pdf(rate_file, m, w)[0]#rate(rate_file, m, w)
        pos = poisson.pmf(x, mu=r)
        #print(f"Rate (mu): {r}")
        # pos = poisson.pmf(x, mu=r) * get_from_pdf(pdf_file, m, w)[0]
        # if is_point_inside_strip(m, w, year, tolerance=0.1)=='True':
        #     pos = poisson.pmf(x, mu=r)
        # else:
        #     pos=0
        # #print(get_from_pdf(pdf_file, m, w)[0])
        return pos
    posterior_data = {}
    
    
    
    for i in range(1, 32):
        print(i)
        # Define the domain for Mosquito (M) and Weather (W) axes
        maxm = get_from_pdf(pdf_file, m, w)[1]
        maxw = get_from_pdf(pdf_file, m, w)[2]
        M = np.linspace(0.5e12,7.3e12, 200)
        W = np.linspace(1800,5100, 200)
        X, Y = np.meshgrid(M, W)
        
        # Compute posterior for each combination of (M, W)
        Z = np.zeros_like(X)
        x_value = i  # The value of observation X = x
        
        for ii in range(X.shape[0]):  # Use 'ii' to avoid overwriting the loop variable 'i'
            for jj in range(X.shape[1]):
                Z[ii, jj] = posterior(X[ii, jj], Y[ii, jj], x_value)
        
        # Normalize the posterior to form a valid PDF
        Z_sum = Z.sum()
        if Z_sum > 0:  # Avoid division by zero
            Z /= Z_sum
        
        # Save posterior data as an array (m, w, density, x_value)
        m_values = X.ravel()
        w_values = Y.ravel()
        density_values = Z.ravel()
        x_values = np.full_like(m_values, fill_value=x_value, dtype=np.float64)
        
        # Combine all values into a single array
        posterior_array = np.column_stack((m_values, w_values, density_values, x_values))
        
        # Save the array to a file and store in the dictionary
        np.save(f'{2007}_posterior_x_{x_value}.npy', posterior_array)  # Save as .npy file
        posterior_data[x_value] = posterior_array  # Store in the dictionary
        
        # Plot the posterior surface with annotations
        plt.figure(figsize=(10, 8))
        
        contour = plt.contourf(X, Y, Z, levels=50, cmap='viridis')
        cb = plt.colorbar(contour)
        # cb.set_label('Rate')
        
        # Plot the posterior as filled contours
        # contour = plt.contourf(X, Y, Z, levels=50, cmap='viridis', alpha=0.7)
        
        # # Add a color bar with a label
        # cbar = plt.colorbar(contour)
        cb.set_label("Posterior Probability", fontsize=15)
        
        # Annotate the maximum posterior point
        max_z = Z.max()
        max_pos = np.unravel_index(Z.argmax(), Z.shape)
        max_m = X[max_pos]
        max_w = Y[max_pos]
        parameter = []
        with open(f'curve{year}.txt', "r") as f:
            for line in f:
                m, w = map(float, line.strip().split())
                parameter.append((m, w))
        param_x = [point[0] for point in parameter]
        param_y = [point[1] for point in parameter]
        #plt.plot(param_x, param_y, 'k--', alpha=0.5, markersize=12, label='Parameters')
        #plt.plot(max_m, max_w, '*k', markersize=8)
        # plt.text(
        #     max_m, max_w, f'{x_value}',
        #     color='white', fontsize=15,
        #     ha='left', va='bottom',
        #     bbox=dict(boxstyle="round,pad=0.3", fc="black", ec="none", alpha=0.7)
        # )
        
        # Title and axis labels
        plt.title(f"Posterior PDF, $\pi(M,W|X={x_value})$", fontsize=15)
        plt.xlabel("f(M)", fontsize=14)
        plt.ylabel("f(W)", fontsize=14)
        
        # Add legend for the maximum posterior point
        plt.legend(fontsize=14)
        
        # Adjust layout and show the plot
        plt.tight_layout()
        # plt.xlim(0.5e12,7.3e12)
        # plt.ylim(1800,5500)
        plt.yticks([])
        plt.savefig(f'posterior-{x_value}.jpg', format='pdf', dpi=1000)
        plt.show()
    
    

 