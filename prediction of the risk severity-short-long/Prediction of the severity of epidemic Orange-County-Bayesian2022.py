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
    
    return result,t


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

def generate_step_function(t, r):
    """
    Given:
      t: array-like of times in days (assumed increasing, length N)
      r: array-like of cumulative or "rate" values at those times (same length),
         assumed to be non-decreasing or mostly so.
    
    This function:
      1) Computes the daily step function (0..max_day),
         where step_values[d] = k means the function has
         crossed integer values up to k by day d.
      2) Aggregates those daily values into 52 weeks for a 1-year period.
         i.e., output[w] is the step value at the end of week w 
         (day 7*w+6), unless that goes beyond the data.

    Returns:
      A list of length 52, weekly_values[w], for w=0..51.
    """

    # 1) Create the daily step function (same code as before)
    t = np.asarray(t, dtype=float)
    r = np.asarray(r, dtype=float)

    max_r = np.floor(r.max())
    if max_r < 1:
        # No crossing of integer 1 => all zeros
        max_day = int(np.floor(t[-1]))
        daily_values = [0]*(max_day + 1)
    else:
        max_r = int(max_r)
        crossing_times = {}

        i = 0
        N = len(t)
        # find crossing time for each integer n
        for n in range(1, max_r+1):
            while i < N-1 and r[i+1] < n:
                i += 1
            if i >= N-1:
                break
            r_i, r_next = r[i], r[i+1]
            t_i, t_next = t[i], t[i+1]
            if r_next == r_i:
                continue
            frac = (n - r_i) / (r_next - r_i)
            crossing_t = t_i + frac*(t_next - t_i)

            # Snap crossing_t to the nearest multiple of 7
            snapped_week = round(crossing_t / 7.0)
            snapped_time = 7.0 * snapped_week
            crossing_times[n] = snapped_time

        max_day = int(np.floor(t[-1]))
        daily_values = [0]*(max_day+1)

        for n in range(1, max_r+1):
            if n not in crossing_times:
                break
            jump_day = int(round(crossing_times[n]))
            if jump_day < 0: 
                continue
            if jump_day > max_day:
                continue
            for d in range(jump_day, max_day+1):
                daily_values[d] = max(daily_values[d], n)

    # 2) Convert daily_values into 52-week output
    #    We'll define "week w" to run from day (7*w) through day (7*w + 6)
    #    and use the value on the last day of the week (7*w + 6),
    #    clamping if beyond max_day.

    weekly_values = []
    for w in range(52):   # 0..51
        day_of_week_end = 7*w + 6  # last day of that week
        if day_of_week_end > len(daily_values)-1:
            day_of_week_end = len(daily_values)-1  # clamp to last day
        weekly_values.append(daily_values[day_of_week_end])

    return weekly_values




#===============================dictionary for saving number of cases and returning year and week of that
#=========================================================================================================
#=========================================================================================================

years = list(range(2006, 2022))  # Years from 2006 to 2020 (excluding 2021)

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

#print("mean",meanC)
#------------------------------------------------------------------------------





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
ending=363#starting_time+8*7


tt=1
print("starting---ending",starting,ending)

years = [year for year in range(2006, 2023) if year != 2010]
print("year",years)
# Dictionary to store cumulative mosquito profiles (M) and weather parameter (W) by year
mosquito_profiles = {}
BR_dic={}

# T = list(np.loadtxt('T-2022.txt'))
# T=T[len(T)-4*365:len(T)-365]
# inuse=T
# model = AutoReg(inuse, lags=365)
# model_fitted = model.fit()
# predictions = model_fitted.predict(start=len(inuse), end=len(inuse)+365)
# TT=list(predictions)
# Dama=TT

# H = list(np.loadtxt('H-2022.txt'))
# H=H[len(H)-4*365:len(H)-365]
# inuseH=H
# modelH = AutoReg(inuseH, lags=365)
# model_fittedH = modelH.fit()
# predictionsH = model_fittedH.predict(start=len(inuseH), end=len(inuseH)+365)
# HH=list(predictionsH)
# Humidity=HH


# P = list(np.loadtxt('P-2022.txt'))
# P=P[len(P)-4*365:len(P)-365]
# inuseP=P
# modelP = AutoReg(inuseP, lags=365)
# model_fittedP = modelP.fit()
# predictionsP = model_fittedP.predict(start=len(inuseP), end=len(inuseP)+365)
# PP=list(predictionsP)
# Precipitation=PP


for year in years: 
    if year==2022:
        d=ending+2
    else:
        d=365
    
    if year==2022:
        dama=get_temperature_for_year(file_path, year)[1][0:d]
        C=meanC[0:d]
        Casetime = list(np.loadtxt(f"./orangecases{year}.txt"))
    elif year==2021:
        dama=get_temperature_for_year(file_path, year)[1]
        C=meanC
        Casetime = list(np.loadtxt(f"./orangecases{year}.txt"))

    else:
        dama = list(np.loadtxt(f"./T{year}.txt"))
        #C=meanC
        C = list(np.loadtxt(f"./Cf{year}.txt"))
        Casetime = list(np.loadtxt(f"./orangecases{year}.txt"))

    # Placeholder for coefficient matrix
    Coematrix = np.zeros((8, 365))  

    # Fill Coematrix based on time-dependent functions

        
    for time in range(0, d):
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
        time = int(divmod(time, d)[1])
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
    if year==2022:
        t = np.linspace(0, d, d-1)
        x0 = (3000000, 0, 0, 0, 10, 10, 10, 10, 10, 100, 100, 100, 100, 100, 100)
        x = odeint(BMPPT, x0, t)
        incB=1/3
        rC=1/6
        muwnvB=0.9
        muB=(0.15/365)

        # Compute cumulative mosquito values (M)
        R=x[0:d,3]
        I=x[0:d,2]
        M=np.cumsum(x[0:d,6])
        Mos=x[0:d,6]
        BS=x[0:d,11]
        TH=x[0:d,0]+x[0:d,1]+x[0:d,2]+x[0:d,3]
        D=x[0:d,11]+x[0:d,12]+x[0:d,13]+x[0:d,14]
        BetaBM=Coematrix[7,0:d-1]*Coematrix[3,0:d-1]/(D+TH)
        BetaMB=Coematrix[3,0:d-1]/(D+TH)
        R2=(BetaBM*BetaMB*Mos*BS*Coematrix[6,0:d-1]*incB)/((muB+incB)*(rC+muwnvB+muB)*(Coematrix[2,0:d-1])*(Coematrix[6,0:d-1]+Coematrix[2,0:d-1]))
        R0=[np.sqrt(R2[i]) for i in range(0,d-1)]
        basicR0=R0
        
        
        nc= generate_step_function(t, R)
        #print(nc)
        newcases=[]
        for i in range(1,len(nc)):
            newcases.append(nc[i]-nc[i-1])
            


       # print("New cases:",newcases)
        # plt.plot(R)
        # plt.show()
        for cases in newcases:  # Week is 1-based
            if cases > 0:  # Only consider weeks with reported casesr
                #print(newcases.index(cases))
                cases_occurrences[cases].append((newcases.index(cases), year))  # Add (week, year) pair
                #print("AAAAAAAAAAAAAAAAAAAAA",cases_occurrences[cases])
        
    else:     
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
    BR_dic[year]=(np.cumsum((get_temperature_for_year(file_path, year)[0])))
    # if year==2022:

    #     Dama = np.array(Dama)
    #     Humidity = np.array(Humidity)
    #     Precipitation = np.array(Precipitation)
        
    #     # Now you can perform arithmetic operations
    #     BR_dic = np.cumsum(0.5 * (Dama + Humidity + Precipitation))
    #     # plt.plot(BR_dic)
    #     # plt.show()
    #     # print("M",len(M))
    #     # print("M",len(BR_dic))
        
    #     Curve2022 = []
    
    #     for i in range(0, 365):
    #         #print(i)
    #         Curve2022.append((mosquito_profiles[i], BR_dic[i]))
    #     #print(Curve2022)
    
    #     # Save to a text file in two columns format
    #     with open("Curve2022_long.txt", "w") as file:
    #         for item in Curve2022:
    #             file.write(f"{item[0]}\t{item[1]}\n")
        
    #     parameter=np.loadtxt('Curve2022_long.txt')


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

# Separate x and y values for all combined datasets
x1, y1 = zip(*all_coordinates)







def kernel_frequency_surface(coordinates, bw=0.3):
    """
    Create a continuous 2D 'frequency' surface over (m, w) coordinates
    using kernel smoothing, where Z represents the frequency at (m, w).

    Parameters:
        coordinates (list of tuples): List of (m, w) coordinates.
        bw (float): Bandwidth for the Gaussian KDE.

    Returns:
        grid_x, grid_y, grid_z: The frequency-based 2D height surface grid.
    """


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
    np.savez("2d_rate_2022.npz", X=grid_x, Y=grid_y, Z=grid_z)
    print("Scaled 2D rate saved as '2d_rate_2022.npz'")



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
    parameter = np.loadtxt("Curve2022.txt")

    param_x = [point[0] for point in parameter]
    param_y = [point[1] for point in parameter]
    plt.scatter(param_x, param_y, color='black', alpha=0.5, s=12, label='Parameters')


    plt.show()

    return grid_x, grid_y, grid_z



grid_x, grid_y, grid_z = kernel_frequency_surface(all_coordinates, bw=0.1)

#=========================================prior==========================================
#=========================================prior==========================================
#=========================================prior==========================================
#=========================================prior==========================================



parameter=np.loadtxt('Curve2022.txt')[starting:ending]
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
    np.savez("2d_prior_test_2022.npz",
              area=area_of_region,
              X=X,
              Y=Y,
              Z=Z)
    print("Saved '2d_prior_test_2022.npz' with keys: [area, X, Y, Z].")

    return area_of_region, shadow, pdf_function


# def plot_curve_with_fixed_shadow_and_uniform_pdf(coordinates, r=100, grid_res=100, alpha=0.5, sigma=None):
# #def plot_curve_with_fixed_shadow_and_uniform_pdf(coordinates, r=100, grid_res=100):
#     """
#     Plots a curve with a "shadow" buffer at distance r, then sets up a PDF over that region.
#     The PDF is a mixture of a uniform density over the shadow and Gaussian "bumps" (2D normals with small variance)
#     centered at a few selected points along the curve.
    
#     Parameters:
#         coordinates (list of (float, float)): The (x,y) points of the curve.
#         r (float): Buffer radius.
#         grid_res (int): Resolution for sampling the PDF on a grid.
#         alpha (float): Weight for the uniform component (0<alpha<1). The bump component gets weight (1-alpha).
#         sigma (float): Standard deviation for the Gaussian bumps. If None, defaults to r/10.
    
#     Returns:
#         area_of_region (float): Area of the buffer region.
#         shadow_polygon (Shapely Polygon): The buffer region.
#         pdf_function (callable): A function pdf(x, y) returning the density at (x,y).
        
#     Also saves a discrete representation of the PDF in "2d_prior_test_bumps.npz" with keys [area, X, Y, Z].
#     """
#     # 1) Build the curve as a Shapely LineString
#     curve = LineString(coordinates)
    
#     # 2) Create the buffer (shadow) around the curve
#     shadow = curve.buffer(r)
    
#     # 3) Compute the area of the shadow
#     area_of_region = shadow.area
    
#     # Set sigma if not provided (small variance relative to the buffer)
#     if sigma is None:
#         sigma = r / 10.0
    
#     # Select a few bump centers along the curve.
#     # For example, choose 5 evenly spaced points along the coordinate list.
#     n = len(coordinates)
#     num_bumps = 5
#     step = max(1, n // num_bumps)
#     bump_centers = [coordinates[i] for i in range(0, n, step)]
    
#     # Define the PDF as a mixture of uniform and Gaussian bump components.
#     def pdf_function(x, y):
#         pt = Point(x, y)
#         if not shadow.contains(pt):
#             return 0.0
#         # Uniform component: density is alpha/area
#         density = alpha / area_of_region
#         # Bump component: equally distribute (1-alpha) weight among all bump centers.
#         bump_weight = (1 - alpha) / len(bump_centers)
#         for center in bump_centers:
#             dx = x - center[0]
#             dy = y - center[1]
#             # 2D normal density: 1/(2*pi*sigma^2) * exp(-((dx^2+dy^2)/(2*sigma^2)))
#             density += bump_weight * (1.0 / (2 * np.pi * sigma**2)) * np.exp(-((dx*dx + dy*dy) / (2 * sigma**2)))
#         return density
    
#     # --- Plot the region, curve, and bump centers ---
#     x_vals, y_vals = zip(*coordinates)
#     shadow_x, shadow_y = shadow.exterior.xy
    
#     fig, ax = plt.subplots(figsize=(10,8))
#     ax.set_facecolor("#8F689A")  # Dark purple background
#     plt.fill(shadow_x, shadow_y, color='green', alpha=0.6,
#              edgecolor='#4B0082', linewidth=2, label='Shadow region')
#     plt.plot(x_vals, y_vals, color='blue', linewidth=2, label='Curve')
    
#     # Mark bump centers
#     bump_x, bump_y = zip(*bump_centers)
#     plt.scatter(bump_x, bump_y, color='red', s=50, zorder=5, label='Bump centers')
    
#     plt.xlabel("X")
#     plt.ylabel("Y")
#     plt.title(f"Curve with Buffer (r={r}) and Gaussian Bumps (σ={sigma}, α={alpha})")
#     plt.legend()
#     plt.show()
    
#     # --- Sample the PDF on a grid ---
#     minx, miny, maxx, maxy = shadow.bounds
#     xs = np.linspace(minx, maxx, grid_res)
#     ys = np.linspace(miny, maxy, grid_res)
#     X, Y = np.meshgrid(xs, ys)
#     Z = np.zeros_like(X)
#     for i in range(grid_res):
#         for j in range(grid_res):
#             Z[i, j] = pdf_function(X[i, j], Y[i, j])
    
#     # Normalize the discrete PDF so that the total mass is 1.
#     dx = (maxx - minx) / (grid_res - 1)
#     dy = (maxy - miny) / (grid_res - 1)
#     total_mass = np.sum(Z) * dx * dy
#     Z /= total_mass
    
#     np.savez("2d_prior_test_bumps.npz", area=area_of_region, X=X, Y=Y, Z=Z)
#     print("Saved '2d_prior_test_bumps.npz' with keys: [area, X, Y, Z].")
    
#     return area_of_region, shadow, pdf_function



area, shadow_polygon, pdf_func = plot_curve_with_fixed_shadow_and_uniform_pdf(parameter, r=200, grid_res=400)
print(f"Area of region = {area}")

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

from scipy.stats import poisson

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

def posterior_array(M_array, W_array, x_value, prior_file="2d_prior_test_2022.npz", rate_file="2d_rate_2022.npz"):
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
    return pmf_values * prior_values

# ------------------------------------------------------------------
# 3) Main code

posterior_data = {}

for i in range(tt, 32):
    pdf_file = "2d_prior_test_2022.npz"

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
                        prior_file="2d_prior_test_2022.npz",
                        rate_file="2d_rate_2022.npz")

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
    np.save(f'posterior_x_{x_value}_2022.npy', posterior_array_4cols)
    posterior_data[x_value] = posterior_array_4cols

    # 3.4) Plot
    plt.figure(figsize=(10, 8))
    contour = plt.contourf(X, Y, Z, levels=50, cmap='viridis', alpha=0.7)
    cbar = plt.colorbar(contour)
    cbar.set_label("Posterior Probability", fontsize=12)

    # Locate the maximum posterior
    max_pos = np.unravel_index(Z.argmax(), Z.shape)
    max_m = X[max_pos]
    max_w = Y[max_pos]

    plt.plot(max_m, max_w, '*k', markersize=8, label=f'Max Posterior (x={x_value})')
    
    
    
    
    parameter = []
    parameter = np.loadtxt("Curve2022.txt")[starting:ending]

    param_x = [point[0] for point in parameter]
    param_y = [point[1] for point in parameter]
    plt.scatter(param_x, param_y, color='black', alpha=0.5, s=12, label='Parameters')
    plt.text(
        max_m, max_w, f'{x_value}',
        color='white', fontsize=10,
        ha='left', va='bottom',
        bbox=dict(boxstyle="round,pad=0.3", fc="black", ec="none", alpha=0.7)
    )

    plt.title("Posterior Probability Surface with Parameters", fontsize=14)
    plt.xlabel("Mosquito Axis (M)", fontsize=12)
    plt.ylabel("Weather Axis (W)", fontsize=12)
    plt.legend(fontsize=10)
    plt.tight_layout()
    plt.show()





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
        filename = f"posterior_x_{x_val}_2022.npy"
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
file_path = 'Curve2022.txt'
loaded_data = np.loadtxt(file_path, skiprows=1)  # Skip the header row

# Extract mosquito and weather profiles
loaded_mosquito_profiles = loaded_data[:, 0]
loaded_weather_profiles = loaded_data[:, 1]

# print(f"Loaded Mosquito Profiles for {year}: {loaded_mosquito_profiles}")
# print(f"Loaded Weather Profiles for {year}: {loaded_weather_profiles}")
# Example usage
x_values = list(range(tt, 32))  # The range of x_values for which posteriors are saved
# m0 = 5.5e12  # Example m-coordinate
# w0 = 4000    # Example w-coordinate

# max_x_value, max_density = find_highest_posterior(x_values, m0, w0)
# print(f"The highest posterior density at ({m0}, {w0}) is for x_value={max_x_value}, with density={max_density}.")
cases=[]
for i in range(starting,ending):
    m=loaded_mosquito_profiles[i]
    w=loaded_weather_profiles[i]
    #print(m,w)
    #print(find_best_x_for_point(m, w, interpolators_dict)[0])
    cases.append(find_best_x_for_point(m, w, interpolators_dict)[0])


    








# def visualize_risk(predicted_cases, mosquito_profile,year,starting,ending):
#     """
#     Visualize the risk of West Nile disease infection in Orange County with mosquito profiles.
    
#     This function creates a figure with two subplots:
    
#       (1) The top subplot shows:
#           - Predicted daily cases as color-coded bars (with colors ranging from yellow to red)
#           - Weekly data (Time) as black stars (assumed one per week, plotted at x = 7 * week_index)
#           - The mosquito profile as a blue line on a secondary y-axis.
          
#     """
#     # Prepare the x-axis for days in the top plot
#     days = np.arange(starting, ending)  # e.g. 0, 1, 2, ..., 363

#     # Create a colormap and normalization for predicted cases (assumes cases from 1 to 32)
#     cmap = plt.colormaps.get_cmap('YlOrRd')  # Yellow-to-red colormap
#     norm = Normalize(vmin=1, vmax=32)


#     fig = plt.figure(figsize=(13, 12))

#     # If you're only using one subplot, you can remove 'height_ratios':
#     gs = fig.add_gridspec(1, 1)  # 1 row, 1 column
    
#     # Create a single subplot (ax1) in the only available grid cell
#     ax1 = fig.add_subplot(gs[0, 0])

#     # ----- TOP PLOT: Original Visualization -----
#     # 1) Plot predicted daily cases as color-coded bars.
#     max_Time = max(Time) if len(Time) > 0 else 1
#     for day, cases in zip(days, predicted_cases[:ending]):
#         bar_color = cmap(norm(cases))
#         ax1.bar(day, height=max_Time, width=1, color=bar_color,
#                 align='edge', edgecolor='none',
#                 label="Predicted Cases (1-32)" if day == days[0] else "")
    
#     # 2) Add a horizontal colorbar for the predicted cases.
#     sm = ScalarMappable(cmap=cmap, norm=norm)
#     cbar = fig.colorbar(sm, ax=ax1, orientation='horizontal', fraction=0.05, pad=0.2)
#     cbar.set_label("Predicted Cases (1 to 32)", fontsize=12)
    
#     # 3) Plot the weekly observed cases as black stars.
#     weeks = np.arange(len(Time))
#     star_x = weeks * 7  # Each week starts at day = 7 * week_index
#     ax1.plot(star_x, Time, '*k', label='Weekly Observed Cases')
    
#     # 4) Add the mosquito profile on a second y-axis.
#     ax2 = ax1.twinx()
#     ax2.plot(np.arange(0, 364), mosquito_profile, color='blue', linewidth=1.5,
#               label="Mosquito Profile")
#     ax2.set_ylabel("Mosquito Profile", fontsize=12, color='blue')
#     ax2.tick_params(axis='y', labelcolor='blue')
    
#     # 5) Set title, x–axis label, and ticks for the top plot.
#     ax1.set_title(f"Severity of West Nile Disease Infection in Orange County with Mosquito Profile ({year})",
#                   fontsize=16)
#     ax1.set_xlabel("Day of the Year", fontsize=12)
#     ax1.set_xticks(np.arange(starting, ending, step=30))
    
#     # 6) Combine legends from ax1 and ax2.
#     lines_labels_1 = ax1.get_legend_handles_labels()
#     lines_labels_2 = ax2.get_legend_handles_labels()
#     lines = lines_labels_1[0] + lines_labels_2[0]
#     labels = lines_labels_1[1] + lines_labels_2[1]
#     ax1.legend(lines, labels, loc='upper left')
    

    
#     plt.tight_layout()
#     plt.show()



def visualize_risk(predicted_cases, mosquito_profile, Casetime, year, starting_point):
    """
    Visualize the risk of infection for West Nile disease in Orange County with mosquito profiles and weekly human cases.

    Parameters:
        predicted_cases (list or array): A vector of predicted cases for 364 days.
        mosquito_profile (list or array): A vector of mosquito profile values for 364 days.
        Casetime (list or array): A vector of weekly human cases (length = 52 weeks).
        year (int): The year for which the data is being visualized.
        starting_point (int): The day from which predictions start.

    Returns:
        None
    """
    if len(predicted_cases) != 364 or len(mosquito_profile) != 364:
        raise ValueError("Both predicted_cases and mosquito_profile must contain 364 values.")
    if len(Casetime) != 52:
        raise ValueError("Casetime must contain 52 weekly values.")
    
    import matplotlib.patches as mpatches  # Import for legend patches
    days = np.arange(0, 364)
    weeks = np.arange(0, 364, step=7)

    # Create colormap and normalization
    cmap = plt.colormaps.get_cmap('YlOrRd')
    norm = Normalize(vmin=1, vmax=32)  # Adjusted to match colorbar label

    fig, ax1 = plt.subplots(figsize=(15, 13))

    # Plot background bars
    for day in days:
        if day < starting_point-5:
            ax1.bar(day, height=max(Casetime), width=1, color='green', 
                    align='edge', edgecolor='none')
        else:
            ax1.bar(day, height=max(Casetime), width=1, 
                    color=cmap(norm(predicted_cases[day])), 
                    align='edge', edgecolor='none')

    # Enhanced colorbar
    sm = ScalarMappable(cmap=cmap, norm=norm)
    cbar = fig.colorbar(sm, ax=ax1, orientation='horizontal', pad=0.15, aspect=50)
    cbar.set_label("Predicted Cases ", fontsize=20, weight='bold')
    cbar.ax.tick_params(labelsize=20)

    # Mosquito profile plot
    ax2 = ax1.twinx()
    ax2.plot(days, mosquito_profile / max(mosquito_profile), color='blue', 
              linewidth=2, label="Mosquito Profile")
    ax2.set_ylabel("Normalized Mosquito Profile", fontsize=14, color='blue', weight='bold')
    ax2.tick_params(axis='y', labelsize=20, labelcolor='blue')

    # Weekly cases plot
    ax1.bar(weeks, Casetime, width=1, color='black', alpha=0.7, 
            label="Weekly Human Cases")
    ax1.set_ylabel("Weekly Human Cases", fontsize=18, color='black', weight='bold')
    ax1.tick_params(axis='y', labelsize=20, labelcolor='purple')

    # Axis formatting
    ax1.set_xlabel("Day of the Year", fontsize=18, weight='bold')
    ax1.set_xticks(np.arange(0, 365, step=30))
    ax1.tick_params(axis='x', labelsize=12)
    ax1.set_title(f"West Nile Disease Risk Assessment - Orange County ({year})", 
                  fontsize=18, weight='bold', pad=20)

    # Legend enhancements
    green_patch = mpatches.Patch(color='green', label='Pre-prediction Period')
    handles, labels = ax1.get_legend_handles_labels()
    handles += ax2.get_legend_handles_labels()[0] + [green_patch]
    labels += ax2.get_legend_handles_labels()[1] + ['Pre-prediction Period']
    
    fig.legend(handles, labels, loc="upper right", 
                bbox_to_anchor=(0.3, 0.85), fontsize=12, 
                title_fontsize=20)

    plt.tight_layout()
    plt.savefig(f'severity{year}.pdf',dpi=1000)
    plt.show()





def calculate_log_likelihood(A, B):
    """
    Compute the log-likelihood score by evaluating how well the observations (B) fit into
    the expected Poisson distribution defined by parameters (A).
    
    Fixes:
    - If a_val == 0 and b_val == 0, log-likelihood = 0 (since P(X=0) = 1).
    - If b_val is outside the 99% domain, assign a log-likelihood of -10.
    - Prevents log(0) by ensuring probabilities are valid.

    Parameters:
        A (list or array): Poisson parameters for each of the 52 weeks.
        B (list or array): Observed values corresponding to each week.

    Returns:
        total_log_likelihood (float): The sum of all log(P) values.
    """
    log_likelihoods = []

    for week_index, (a_val, b_val) in enumerate(zip(A, B)):
        # Compute 99% quantile for Poisson(a_val)
        quantile_99 = int(poisson.ppf(0.99, a_val))  # Ensure integer for domain limit

        # Handle special case: a_val == 0 and b_val == 0
        if a_val == 0 and b_val == 0:
            log_likelihoods.append(0)  # log(1) = 0
            continue

        # If b_val is outside the domain (greater than 99% quantile), assign log-likelihood = -10
        if b_val > quantile_99:
            log_likelihoods.append(-10)
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
            if b_val < bin_edges[0]:
                bin_index = 0
            else:
                bin_index = len(bin_edges) - 2  # Assign to last bin

        # Compute probability for the bin
        P_bin = poisson.cdf(bin_edges[bin_index + 1], a_val) - poisson.cdf(bin_edges[bin_index], a_val)

        # Avoid log(0) errors but check for perfect probability
        if np.isclose(P_bin, 1.0):
            log_likelihoods.append(0)  # log(1) = 0
        else:
            #print(np.log(max(P_bin, 1e-10)))
            log_likelihoods.append(np.log(max(P_bin, 1e-10)))  # Ensure P_bin is valid

    # Sum all log(P) values
    total_log_likelihood = sum(log_likelihoods)

    return total_log_likelihood









def plot_weekly_data(A, B):
    """
    Plot weekly predicted data (A) vs. weekly observed data (B).
    
    Parameters:
        A (list): Weekly predicted cases.
        B (list): Weekly real reported cases.
    """
    weeks = np.arange(len(A))  # Week indices: 0, 1, 2, ..., 51
    x_vals = weeks * 7  # Convert weeks into day positions (7-day intervals)

    plt.figure(figsize=(13, 6))

    # Plot predicted cases (A) as a green line with markers
    plt.plot(x_vals, A, 'g-o', label='Predicted Cases')

    # Plot observed cases (B) as a red line with markers
    plt.plot(x_vals, B, 'r-o', label='Reported Cases')

    # Title and labels
    plt.title("Weekly Predicted vs. Reported Cases")
    plt.xlabel("Day of the Year (7-day intervals)")
    plt.ylabel("Number of Cases")

    # Set x-axis labels as "Week 1", "Week 2", etc.
    plt.xticks(x_vals, [f"Week {w+1}" for w in weeks], rotation=45)

    plt.grid(True)
    plt.legend()
    plt.tight_layout()
    plt.show()

# Generate weekly data





MM=list(np.loadtxt("./MMosquito_Profile_2021.txt"))

#visualize_risk(cases, MM, 2022)


predicted_case=[0 for i in range(0,starting+1)]
predicted_case=predicted_case+cases

print(len(cases))
print("Length of predicted_case:", len(predicted_case))
print("Length of MM (mosquito_profile):", len(MM))



B=Time
A1=[0 for i in range(0,int(starting/7))]
A2=[cases[i] for i in range(len(cases)) if i % 7 == 0]  
A=A1+A2
#visualize_risk(cases, MM, year, starting, ending)
visualize_risk(predicted_case, MM, B,year, starting)




print(calculate_log_likelihood(A[int(starting/7):], B[int(starting/7):]))
# Plot the data
plot_weekly_data(A, B)

AA=[0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 1, 1, 2, 2, 2, 3, 3, 3, 3, 2, 2, 2, 1, 1, 1, 0, 0, 0, 0]#AA=[0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 1, 1, 2, 2, 3, 3, 3, 3, 3, 3, 3, 2, 2, 1, 1, 1, 0, 0, 0]
print(calculate_log_likelihood(AA[int(starting/7):], B[int(starting/7):]))

x=calculate_log_likelihood(A[int(starting/7):], B[int(starting/7):])
y=calculate_log_likelihood(AA[int(starting/7):], B[int(starting/7):])

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
    plt.figure(figsize=(10, 6))
    
    # Plot real cases as black bars (thicker)
    plt.bar(weeks, B, color='black', alpha=0.7, label='Reported Values', width=0.6)

    # Plot model results as lines with scores in the legend
    plt.plot(weeks, A, marker='o', linestyle='-', label=f'Probabilistic Model Result (Score = {np.round(x,decimals=2)})', linewidth=2)
    plt.plot(weeks, AA, marker='s', linestyle='--', label=f'Negative Binomial Model (Score = {np.round(y,decimals=2)})', linewidth=2)

    # Labels and title
    plt.xlabel('Week')
    plt.ylabel('Values')
    plt.title(f'Comparison of Model Results and Reported Cases - {year}')
    plt.legend()
    
    # Save as PDF
    filename = f"model_comparison_{year}.pdf"
    plt.savefig(filename, format="pdf")

    # Show the plot
    plt.show()

    print(f"Plot saved as {filename}")
year=2022


plot_model_comparison(A, AA, B, year,x,y)