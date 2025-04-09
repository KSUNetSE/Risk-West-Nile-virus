import numpy as np
import matplotlib.pyplot as plt
from scipy.integrate import odeint
from matplotlib import pyplot as plt
from scipy.stats import gamma
from scipy.stats import gaussian_kde
import pandas as pd
from scipy.optimize import minimize
from scipy.stats import gamma, norm, multivariate_normal, chi2, kstest



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


#years = [year for year in range(2006, 2021) if year != 2010]
# Dictionary to store cumulative mosquito profiles (M) by year
mosquito_profiles = {}
Weather={}

        
 
#tem2021=get_temperature_for_year(file_path, 2021)
    
years=[2020]
for year in years:
    print("year",year)
    if year==2021:
        dama=T2021
        C = Cf2021
    else:
        dama = list(np.loadtxt(f"./T{year}.txt"))
        #H = list(np.loadtxt(f"./H{year}.txt"))
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
    Weather[year]=np.cumsum(get_temperature_for_year(file_path, year))[0:364]
    temperatures = get_temperature_for_year(file_path, year)
    if  year == 2021:
        # Extract data for the year
        mosquito_profile_2021 = mosquito_profiles[year]
        weather_profile_2021 = Weather[year]
        
        # Combine into a single array for saving
        data_to_save = np.column_stack((mosquito_profile_2021, weather_profile_2021))
        
        # Save to a .txt file
        np.savetxt(f"MMosquito_Profile_{year}.txt", Mos, header="Mosquito_Profile", fmt='%.6f')

        
        np.savetxt(f"Weather_Mosquito_Profiles_{year}.txt", data_to_save, header="Mosquito_Profiles\tWeather_Profiles", fmt='%.6f')
        print(f"Weather and mosquito profiles for {year} saved to 'Weather_Mosquito_Profiles_{year}.txt'")

            
    
#     plt.plot(mosquito_profiles[year],Weather[year],'r.')
# plt.show()
    



# Generate sample data
    daybegin=0
    dayend=364
    parameter = []
    curve=[]
    for i in range(daybegin, dayend):
        point_mean = Weather[year][i]  # Replace with your Weather[year][i] array
        samples = np.random.normal(point_mean, scale=1, size=30)
        curve.append((mosquito_profiles[year][i], Weather[year][i]))
        
        
        
        for y in samples:
            parameter.append((mosquito_profiles[year][i], y))  # Replace mosquito_profiles[year][i]
            
    
    with open(f'Curve{year}.txt', "w") as f:
        for param in curve:
            f.write(f"{param[0]}\t{param[1]}\n")

    print("Parameters saved to",f'Curve{year}.txt')
    
    # Extract x and y values from parameter
    x = [point[0] for point in parameter]
    y = [point[1] for point in parameter]
    
    # Estimate density using KDE
    data = np.vstack([x, y])
    kde = gaussian_kde(data)
    density = kde(data)
    
    # Plot the data with density-based coloring
    #plt.figure(figsize=(10, 8))
    plt.scatter(x, y, c=density, cmap='viridis', s=10, alpha=0.7)
    # plt.colorbar(label='Density')
    # plt.title('Density-Based Coloring of Samples')
    # plt.xlabel('Mosquito Profiles (X-axis)')
    # plt.ylabel('Weather Profiles (Y-axis)')
    # plt.tight_layout()
    
    
    
    # for i in range(daybegin, dayend):
    #     point_mean = Weather[year][i]
    #     samples = np.random.normal(point_mean, scale=100, size=50)
    
    #     for y in samples:
    #         parameter.append((mosquito_profiles[year][i], y))

    # # Extract x and y values from the combined parameter list
    # x = [point[0] for point in parameter]
    # y = [point[1] for point in parameter]
    
    # # Estimate density using KDE for all data
    # data = np.vstack([x, y])
    # kde = gaussian_kde(data)
    # density = kde(data)
    
    # # Plot the combined data with density-based coloring on one page

    # plt.scatter(x, y, c=density, cmap='viridis', s=10, alpha=0.7)

    
# plt.colorbar(label='Density')
# plt.figure(figsize=(12, 10))
# plt.title('Density-Based Coloring of All Samples (Single Plot)')
# plt.xlabel('Mosquito Profiles (X-axis)')
# plt.ylabel('Weather Profiles (Y-axis)')
# #plt.grid(alpha=0.3)
# plt.tight_layout()
plt.tick_params(axis='y', labelleft=False)
plt.show()







        





import numpy as np
import matplotlib.pyplot as plt
from scipy.interpolate import interp1d
from scipy.stats import norm

def generate_2d_normal_pdf(coordinates, year, num_samples=10000, sigma=1.0):
    """
    Generate a 2D normal distribution around a curve defined by a list of coordinates.

    Parameters:
        coordinates (list of tuples): List of (x, y) coordinates representing the curve.
        year (int or str): Year or identifier for saving the file.
        num_samples (int): Number of random samples to generate.
        sigma (float): Standard deviation of the normal distribution (controls spread).

    Returns:
        X, Y, Z: Arrays representing the grid and the normal PDF.
    """
    # Extract x and y values from the coordinates
    x_values = np.array([point[0] for point in coordinates])
    y_values = np.array([point[1] for point in coordinates])

    # Interpolate the curve to approximate y as a function of x
    curve_interp = interp1d(x_values, y_values, kind='cubic', fill_value="extrapolate")

    # Generate a grid for the PDF
    x_grid = np.linspace(min(x_values), max(x_values), 100)
    y_grid = np.linspace(min(y_values) - 3 * sigma, max(y_values) + 3 * sigma, 100)
    X, Y = np.meshgrid(x_grid, y_grid)

    # Compute the normal PDF centered around the curve
    Z = np.exp(-0.5 * ((Y - curve_interp(X)) / sigma) ** 2) / (sigma * np.sqrt(2 * np.pi))

    # Normalize the PDF so that it integrates to 1 over the grid
    Z /= np.sum(Z) * (x_grid[1] - x_grid[0]) * (y_grid[1] - y_grid[0])

    # Plot the results
    plt.figure(figsize=(10, 8))
    plt.plot(x_values, y_values, '-', label='C(M,W)', color='blue')
    plt.contourf(X, Y, Z, levels=50, cmap='viridis', alpha=0.6, label='Normal PDF')
    #plt.colorbar(label='Density')
    plt.xlabel('f(M)',fontsize=14)
    plt.ylabel('f(W)',fontsize=14)
    plt.yticks([])
    plt.xticks([])
    plt.title(f'2D Normal Prior based on (M,W)',fontsize=15)
    plt.legend(loc='upper left')
    plt.savefig('prior_Normal.pdf', format='pdf', bbox_inches='tight', dpi=1000)

    plt.show()

    # Save the data as a .npz file
    #np.savez(f'2d_prior_pdf-{year}.npz', X=X, Y=Y, Z=Z)

    return X, Y, Z

# year=2007
# generate_2d_normal_pdf(curve, year, num_samples=10000, sigma=100)
    


import numpy as np
import matplotlib.pyplot as plt
from scipy.interpolate import interp1d



def generate_2d_uniform_pdf(coordinates, year, strip_width=30):
    """
    Generate a 2D uniform distribution around a curve and save only the (x, y) points in the strip.

    Parameters:
        coordinates (list of tuples): List of (x, y) coordinates representing the curve.
        year (int or str): Year or identifier for saving the file.
        strip_width (float): Width of the strip around the curve.

    Returns:
        strip_x, strip_y: Arrays of (x, y) points inside the uniform strip.
    """
    # Extract x and y values from the coordinates
    x_values = np.array([point[0] for point in coordinates])
    y_values = np.array([point[1] for point in coordinates])

    # Interpolate the curve to approximate y as a function of x
    curve_interp = interp1d(x_values, y_values, kind='cubic', fill_value="extrapolate")

    # Generate a grid for the PDF
    x_grid = np.linspace(min(x_values), max(x_values), 100)
    y_grid = np.linspace(min(y_values) - strip_width, max(y_values) + strip_width, 100)
    X, Y = np.meshgrid(x_grid, y_grid)

    # Compute the uniform PDF within the strip around the curve
    strip_mask = (Y >= curve_interp(X) - strip_width) & (Y <= curve_interp(X) + strip_width)

    # Extract (x, y) coordinates where the mask is True
    strip_x = X[strip_mask]
    strip_y = Y[strip_mask]

    # Save only the (x, y) points in the uniform strip
    np.savez(f'2d_prior_pdf-{year}.npz', strip_x=strip_x, strip_y=strip_y)

    # Plot the results
    plt.figure(figsize=(8, 6))
    plt.plot(x_values, y_values, '-', label='Input Curve', color='blue')
    plt.scatter(strip_x, strip_y, s=1, color='green', alpha=0.6, label='Uniform Strip Area')
    plt.xlabel('x')
    plt.ylabel('y')
    plt.title(f'2D Uniform Strip around Input Curve ({year})')
    plt.legend()
    plt.grid(True)
    plt.show()

    return strip_x, strip_y

import numpy as np
import matplotlib.pyplot as plt
import seaborn as sns
from scipy.interpolate import interp1d

def generate_2d_uniform_pdf_improved(coordinates, year, strip_width=1000):
    """
    Generate a 2D uniform distribution around a curve and visualize density more effectively.

    Parameters:
        coordinates (list of tuples): List of (x, y) coordinates representing the curve.
        year (int or str): Year or identifier for saving the file.
        strip_width (float): Width of the strip around the curve.

    Returns:
        strip_x, strip_y: Arrays of (x, y) points inside the uniform strip.
    """
    # Extract x and y values from the coordinates
    x_values = np.array([point[0] for point in coordinates])
    y_values = np.array([point[1] for point in coordinates])

    # Interpolate the curve to approximate y as a function of x
    curve_interp = interp1d(x_values, y_values, kind='cubic', fill_value="extrapolate")

    # Generate a grid for the PDF
    x_grid = np.linspace(min(x_values), max(x_values), 100)
    y_grid = np.linspace(min(y_values) - strip_width, max(y_values) + strip_width, 100)
    X, Y = np.meshgrid(x_grid, y_grid)

    # Compute the uniform PDF within the strip around the curve
    strip_mask = (Y >= curve_interp(X) - strip_width) & (Y <= curve_interp(X) + strip_width)

    # Extract (x, y) coordinates where the mask is True
    strip_x = X[strip_mask]
    strip_y = Y[strip_mask]

    # Save only the (x, y) points in the uniform strip
    np.savez(f'2d_prior_pdf-{year}.npz', strip_x=strip_x, strip_y=strip_y)

    # Create a density plot
    plt.figure(figsize=(10, 8))
    plt.plot(x_values, y_values, '-', label='Input Curve', color='blue', linewidth=2)
    sns.kdeplot(x=strip_x, y=strip_y, cmap="viridis", fill=True, alpha=0.7, levels=50)
    
    plt.xlabel('x')
    plt.ylabel('y')
    plt.title(f'2D Density of Uniform Strip around Input Curve ({year})')
    plt.legend()
    plt.xlim(0,6.5e12)
    plt.show()

    return strip_x, strip_y




from shapely.geometry import LineString
import geopandas as gpd
def plot_curve_with_fixed_shadow(coordinates, r=100):
    """
    Plots a curve with a shadow effect at a fixed distance and highlights the zero-density region.

    Parameters:
        coordinates (list of tuples): List of (x, y) coordinates representing the curve.
        r (float): Radius of the shadow effect.
    """
    # Convert the coordinates into a LineString
    curve = LineString(coordinates)

    # Create a shadow (buffer around the curve)
    shadow = curve.buffer(r)

    # Extract x and y values of the curve
    x_values, y_values = zip(*coordinates)

    # Extract the exterior ring of the buffer polygon for plotting
    shadow_x, shadow_y = shadow.exterior.xy

    # Plot
    fig, ax = plt.subplots(figsize=(10,8))

    # Fill the entire background with dark purple to represent 0 density
    ax.set_facecolor("#8F689A")  # Dark purple

    # Plot shadow region with green fill and dark purple border
    plt.fill(shadow_x, shadow_y, color='green', alpha=0.6, edgecolor='#4B0082', linewidth=2, label='Uniform density')

    # Plot original curve
    plt.plot(x_values, y_values, color='blue', linewidth=2, label='(M,W)')



    plt.xlabel("f(M)",fontsize=14)
    plt.yticks([])
    plt.ylabel("f(W)", fontsize=14)
    plt.xticks([])
    plt.title("Uniform Prior based on C(M, W)",fontsize=15)
    plt.legend()
    plt.savefig('prior_uniform22.pdf', format='pdf', bbox_inches='tight', dpi=1000)
    

    plt.show()
    






coordinates =curve[100:365]
# Plot with the corrected shadow effect
plot_curve_with_fixed_shadow(coordinates, r=300)



generate_2d_normal_pdf(coordinates, year, num_samples=10000, sigma=100.0)





    


def sample_prior(parameters, bandwidth=0.5, num_samples=1):
    """
    Fit a 2D PDF using KDE for a given list of (m, w) coordinates,
    and then sample from this estimated PDF, enforcing constraints on m and w.

    Parameters:
        parameters (list of tuples): List of (m, w) coordinates.
        bandwidth (float): Bandwidth for the KDE.
        num_samples (int): Number of samples to draw from the KDE.

    Returns:
        X, Y, Z: Grid coordinates and the fitted PDF values on the grid.
        samples: Array of shape (num_samples, 2) containing (m, w) samples.
    """
    # Convert the list of (m, w) tuples into two separate arrays
    m_values, w_values = zip(*parameters)
    m_values = np.array(m_values)
    w_values = np.array(w_values)

    # Enforce bounds on input data
    m_values = np.clip(m_values, 0, 8.1e12)
    w_values = np.clip(w_values, 0, 9000)

    # Stack m and w values for KDE
    data = np.vstack([m_values, w_values])  # Shape: (2, num_points)
    kde = gaussian_kde(data, bw_method=bandwidth)

    # Sample from the KDE
    sampled_data = kde.resample(num_samples)  # shape: (2, num_samples)
    # Reformat to (num_samples, 2) and enforce constraints:
    samples = np.column_stack((sampled_data[0], sampled_data[1]))
    samples[:, 0] = np.clip(samples[:, 0], 0, 8.1e12)  # Clamp m
    samples[:, 1] = np.clip(samples[:, 1], 0, 9000)    # Clamp w

    # Create a grid for visualization only (not strictly required for sampling)
    m_grid = np.linspace(0, 7.1e12, 200)
    w_grid = np.linspace(0, 6000, 200)
    X, Y = np.meshgrid(m_grid, w_grid)
    grid_positions = np.vstack([X.ravel(), Y.ravel()])

    # Evaluate the KDE on the grid (for plotting)
    Z = kde(grid_positions).reshape(X.shape)

    # Optionally, save the 2D PDF data
    np.savez("2d_prior_pdf.npz", X=X, Y=Y, Z=Z)
    #np.savez("2d_pdf.npz", X=X, Y=Y, Z=Z)

    # Plot the 2D PDF
    plt.figure(figsize=(10, 8))
    plt.contourf(X, Y, Z, levels=50, cmap='viridis')
    plt.colorbar(label='PDF Density')
    plt.title('2D PDF Fitted with KDE')
    plt.xlabel('Mosquito Axis (m)')
    plt.ylabel('Weather Axis (w)')

    # Optionally, plot the sampled points on top
    plt.scatter(samples[:, 0], samples[:, 1], color='red', alpha=0.3, s=5, label='Samples')
    plt.legend()
    plt.tight_layout()
    plt.show()

    return X, Y, Z, samples

  



import matplotlib.pyplot as plt
from shapely.geometry import LineString, Point
import numpy as np


def plot_curve_with_fixed_shadow_and_uniform_pdf(coordinates, r=100):
    """
    Plots a curve with a "shadow" (buffer) at distance r, 
    sets up a uniform PDF on that region, and returns:
      1. area_of_region (float),
      2. region_polygon (Shapely Polygon),
      3. pdf_function (callable), i.e. pdf(x,y).

    Parameters:
        coordinates (list of (float,float)): The (x,y) points of the curve.
        r (float): Radius of the buffer around the curve.
    """

    # 1) Build the curve as a Shapely LineString
    curve = LineString(coordinates)

    # 2) Create the shadow region (buffer around the curve)
    shadow = curve.buffer(r)

    # 3) Compute the area
    area_of_region = shadow.area

    # 4) Define a uniform PDF over that region
    #    pdf(x,y) = 1/area if inside shadow, else 0
    def pdf_function(x, y):
        pt = Point(x, y)
        if shadow.contains(pt):
            return 1.0 / area_of_region
        else:
            return 0.0

    # --- (Optional) Plotting code below ---
    x_values, y_values = zip(*coordinates)
    shadow_x, shadow_y = shadow.exterior.xy

    fig, ax = plt.subplots(figsize=(10,8))
    ax.set_facecolor("#8F689A")  # dark purple

    plt.fill(shadow_x, shadow_y, color='green', alpha=0.6,
             edgecolor='#4B0082', linewidth=2, label='Uniform region')

    plt.plot(x_values, y_values, color='blue', linewidth=2, label='Curve (M,W)')

    plt.xlabel("M")
    plt.ylabel("W")
    plt.title(f"Uniform PDF Based on C(M,W)",fontsize=14)
    plt.legend()
    plt.savefig('prior_uniform1.pdf', format='pdf', bbox_inches='tight', dpi=1000)
    plt.show()
    # --- End Plot ---

    # -------------------------------------------------------------
    # NEW: Save the area to a NumPy file named "2d_prior_test.npy"
    # -------------------------------------------------------------


    # Return the three items
    px, py = shadow.exterior.xy
    np.savez("2d_prior_test.npz", area=area_of_region, poly_x=px, poly_y=py)
    return area_of_region, shadow, pdf_function

coords =coordinates
plot_curve_with_fixed_shadow_and_uniform_pdf(coords, r=300)



