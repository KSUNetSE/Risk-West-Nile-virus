import numpy as np
import matplotlib.pyplot as plt
from scipy.integrate import odeint
from matplotlib import pyplot as plt
from scipy.stats import gamma
from scipy.stats import gaussian_kde
import pandas as pd
from scipy.optimize import minimize
from scipy.stats import gamma, norm, multivariate_normal, chi2, kstest


#from scipy.stats import norm, multivariate_normal, kstest
#from copulas.visualization import scatter_2d
#from copulas.multivariate import GaussianMultivariate

# Initialize an array to store the total cases for each week across all years
weeks_total_cases = np.zeros(52)  # Assuming 52 weeks in a year

# List of years to process
#years = list(range(2006, 2021))  # Years from 2006 to 2020 (excluding 2021)
years = [year for year in range(2006, 2021) if year != 2010]
specific_year = 2014  # Year to compare
specific_year_cases = None  # Placeholder for the specific year's data

# Read case data for each year
for year in years:
    try:
        case_data = np.loadtxt(f"./orangecases{year}.txt")
        if year == specific_year:
            specific_year_cases = case_data  # Save the specific year's data
        for week, cases in enumerate(case_data):
            weeks_total_cases[week] += cases  # Accumulate cases per week
    except FileNotFoundError:
        print(f"Data for {year} not found. Skipping this year.")

# Plot histogram for cumulative weekly cases
plt.figure(figsize=(12, 8))
weeks = np.arange(1, 53)  # Week numbers from 1 to 52
plt.bar(weeks, weeks_total_cases, color='blue', alpha=0.7, edgecolor='black', label="Cumulative Cases (2006-2020)")

# Overlay specific year's data
if specific_year_cases is not None:
    plt.plot(weeks, specific_year_cases, '-ro', label=f"Cases in {specific_year}", linewidth=2)

# Add labels and title
plt.xlabel("Week Number", fontsize=14)
plt.ylabel("Total Cases (Frequency)", fontsize=14)
plt.title(f"Total Cases per Week (2006-2020) with {specific_year} Overlay", fontsize=16)

# Add legend
plt.legend()

# Improve x-axis visibility
plt.xticks(range(1, 53, 2), rotation=45)  # Show every second week for better readability

# Display the plot
plt.tight_layout()
plt.show()



# List of years to process
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



# Dictionary to store cumulative mosquito profiles (M) by year
mosquito_profiles = {}
BR_dic={}

# List of years to process
years = [year for year in range(2006, 2021) if year != 2010]


for year in years:
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

    
    BR_dic[year]=(np.cumsum((get_temperature_for_year(file_path, year))))



    
    temperatures = get_temperature_for_year(file_path, year)
    
#--------------------------------------------------------------------profile weather plotting    
    # fig, ax1 = plt.subplots(figsize=(10, 6))
    
    # # Plot Casetime on the primary y-axis
    # ax1.plot([i * 7 for i in range(1, 53)], Casetime, '--r', label='Casetime')
    # ax1.set_xlim(100, 365)
    # ax1.set_xlabel('Days')
    # ax1.set_ylabel('Casetime', color='red')
    # ax1.tick_params(axis='y', labelcolor='red')
    
    # # Create a secondary y-axis for temperature
    # ax2 = ax1.twinx()
    
    # # Plot the moving average of temperatures
    # moving_avg_temp = moving_average(temperatures, window_size=1)
    # ax2.plot(range(len(moving_avg_temp)), moving_avg_temp, '--k', label='Temperature (Moving Average)')
    # ax2.set_ylabel('Temperature', color='black')
    # ax2.tick_params(axis='y', labelcolor='black')
    # ax2.set_ylim(0, 40)
    # # Add a title and legend
    # plt.title(f'{year}')
    # fig.legend(loc="upper right", bbox_to_anchor=(0.85, 0.85))
    # fig.tight_layout()
    
    # # Show the plot
    # plt.show()




#--------------------------------------------------------------------------------

data=[]
dataBR=[]

data1=[]
dataBR1=[]

data2=[]
dataBR2=[]


#--------------------------------------------------------------------------
Allmos=[]
Allbir=[]
for i in range(1, 32):
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
                mosquito_value = mosquito_profiles[year][7 * (week-1)]  # Week is 1-based, adjust for 0-indexed array
                BR_value=BR_dic[year][7 * (week-1)]
            Allmos.append(mosquito_value)
            Allbir.append(BR_value)

#------------------------------------------------------------------------------

border=13

A=[]
B=[]
for i in range(1,32):
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
            
        print("moss",mosquito_values_list)
        
        # Plotting the mean mosquito value for this num_cases
        if mosquito_values_list:  # Check if the list is not empty
            if i<border:
                #plt.plot(i, np.mean(mosquito_values_list), 'r*')  # Plot the mean for each num_cases
                plt.plot([i for j in range(0,len(mosquito_values_list))], mosquito_values_list, '.g')
                #plt.plot([i for j in range(0,len(Bird_value_list))], Bird_value_list, '.r')
            else:
                #plt.plot(i, np.mean(mosquito_values_list), 'k*')  # Plot the mean for each num_cases
                plt.plot([i for j in range(0,len(mosquito_values_list))], mosquito_values_list, '.b')
                #plt.plot([i for j in range(0,len(Bird_value_list))], Bird_value_list, '.k')
                
            # plt.plot(i, np.mean(dataBR1), 'r*')  # Plot the mean for each num_cases
            # plt.plot([i for j in range(0,len(dataBR1))], dataBR1, 'bo')


    else:
        print(f"No data found for {num_cases} cases. Skipping.")

#Display the plot
plt.xlabel("Number of Cases")
plt.ylabel("Mosquito Value")
plt.title("Mosquito Values for Different Cases")
plt.grid(True)
plt.show()

plt.plot(data1,dataBR1,'g*')
plt.plot(data2,dataBR2,'r.')
# plt.ylim(0,1)
plt.show()


#--------------------------------------------------------------------------------------






# Combine datasets into coordinate pairs
coordinates1 = list(zip(data1, dataBR1))
coordinates2 = list(zip(data2, dataBR2))

# Combine both datasets into one
all_coordinates = coordinates1 + coordinates2

# Separate x and y values for all combined datasets
x1, y1 = zip(*all_coordinates)

# Create 2D histogram
plt.figure(figsize=(10, 8))
hist, xedges, yedges, im = plt.hist2d(x1, y1, bins=20, cmap='Reds', alpha=0.6, label='Dataset 1')

# Add colorbar for the histogram
plt.colorbar(label='Dataset 1 Density (Reds)', orientation='vertical')

# Fit a PDF using Gaussian KDE
data = np.vstack([x1, y1])  # Combine x and y for KDE
kde = gaussian_kde(data)    # Fit KDE to the data

# Evaluate the KDE on a grid
x_grid = np.linspace(xedges[0], xedges[-1], 200)
y_grid = np.linspace(yedges[0], yedges[-1], 200)
X, Y = np.meshgrid(x_grid, y_grid)
positions = np.vstack([X.ravel(), Y.ravel()])
Z = kde(positions).reshape(X.shape)



# Save the PDF data (X, Y, Z) to a file
# np.savez("2d_pdf.npz", X=X, Y=Y, Z=Z)
# print("PDF saved as '2d_pdf.npz'")

# Overlay the PDF as filled contours (colored by density)
contour = plt.contourf(X, Y, Z, levels=50, cmap='viridis', alpha=0.7)

# Add a colorbar for the contours
plt.colorbar(contour, label='Fitted PDF Density')

# Add grid borders to match histogram bins
ax = plt.gca()  # Get the current axis
for x in xedges:
    ax.axvline(x, color='gray', linestyle='--', linewidth=0.5)  # Vertical gridlines
for y in yedges:
    ax.axhline(y, color='gray', linestyle='--', linewidth=0.5)  # Horizontal gridlines

#Scatter plot of points
plt.scatter(x1, y1, color='black', alpha=0.5, s=10, label='Data Points')

#Customize the plot
plt.xlabel('X')
plt.ylabel('Y')
plt.title('2D Histogram with Colored PDF Contours')
plt.legend()
plt.show()


# Combine datasets into coordinate pairs
coordinates1 = list(zip(data1, dataBR1))
coordinates2 = list(zip(data2, dataBR2))

# Combine both datasets into one
all_coordinates = coordinates1 + coordinates2









from collections import Counter
from scipy.stats import gaussian_kde
import numpy as np
import matplotlib.pyplot as plt
from collections import Counter
from scipy.stats import gaussian_kde





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

    # Force the domain to start at 0
    x_min = 0.0  
    x_max = 9.1e12
    y_min = 0
    y_max = 6000

    # Create a grid for evaluation
    n_points = 400
    x_grid = np.linspace(x_min, x_max, n_points)
    y_grid = np.linspace(y_min, y_max, n_points)
    grid_x, grid_y = np.meshgrid(x_grid, y_grid)
    grid_positions = np.vstack([grid_x.ravel(), grid_y.ravel()])

    # Evaluate the smoothed frequency surface
    grid_z = kde(grid_positions).reshape(grid_x.shape)

    # Find the maximum density and scale the output
    maxden = grid_z.max()
    print("Maximum density (maxden):", maxden)
    scaling_factor = 31 / maxden
    grid_z *= scaling_factor  # Scale the height surface
    np.savez("2d_rate.npz", X=grid_x, Y=grid_y, Z=grid_z)
    print("Scaled 2D rate saved as '2d_rate.npz'")

    # Plot the frequency surface using imshow (no contours)
    plt.figure(figsize=(8, 5))
    plt.imshow(grid_z, extent=(x_min, x_max, y_min, y_max), origin='lower',
               cmap='viridis', aspect='auto')
    plt.colorbar(label='Rate')
    contour_lines = plt.contour(grid_x, grid_y, grid_z, levels=15, colors='black', linewidths=0.3)
    plt.clabel(contour_lines, inline=True, fontsize=8)  # Add labels to the contour lines

    # Overlay scatter plot of original points
    #plt.scatter(x_vals, y_vals, color='red', alpha=0.7, s=10, label='Observations')

    plt.title('Rate of the Poisson Process', fontsize=15)
    plt.xlabel('(M)', fontsize=14)
    plt.ylabel('(W)', fontsize=14)
    plt.tick_params(axis='y', labelleft=False)
    plt.legend()
    plt.tight_layout()
    plt.xlim(1e12, 7e12)
    plt.ylim(2000, 5000)
    plt.savefig('rate_plot.pdf', format="pdf", dpi=1000, bbox_inches='tight')
    plt.show()

    return grid_x, grid_y, grid_z




grid_x, grid_y, grid_z = kernel_frequency_surface(all_coordinates, bw=0.4)







