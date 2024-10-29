# -*- coding: utf-8 -*-
"""
Created on Mon May 20 12:35:23 2024

@author: shosseini
"""

import numpy as np
from scipy.integrate import odeint
from matplotlib import pyplot as plt
from scipy.optimize import curve_fit
import pandas as pd

def get_temperature_for_year(file_path, year):
    # Read the CSV file
    df = pd.read_csv(file_path)
    df_year = df[df['Year'] == year]
    # Convert the temperature column to a list of real values (floats)
    t = df_year['T'].tolist()
    h = df_year['H'].tolist()
    p = df_year['P'].tolist()
    return t,h,p
file_path = 'losAngelesinfo.csv'  # Replace with your actual file path



Cf2006 = list(np.loadtxt("./C2006.txt"))
Cf2007 = list(np.loadtxt("./C2007.txt"))
Cf2008 = list(np.loadtxt("./C2008.txt"))
Cf2009 = list(np.loadtxt("./C2009.txt"))
Cf2010 = list(np.loadtxt("./C2010.txt"))
Cf2011 = list(np.loadtxt("./C2011.txt"))
Cf2012 = list(np.loadtxt("./C2012.txt"))
Cf2013 = list(np.loadtxt("./C2013.txt"))
Cf2014 = list(np.loadtxt("./C2014.txt"))
Cf2015 = list(np.loadtxt("./C2015.txt"))
Cf2016 = list(np.loadtxt("./C2016.txt"))
Cf2017 = list(np.loadtxt("./C2017.txt"))
Cf2018 = list(np.loadtxt("./C2018.txt"))
Cf2019 = list(np.loadtxt("./C2019.txt"))
Cf2020 = list(np.loadtxt("./C2020.txt"))
#--------------------------------------Section of importing data 
#C2000 = list(np.loadtxt("./C2000.txt"))
T2000 = get_temperature_for_year(file_path, 2000)[0]
H2000 = get_temperature_for_year(file_path, 2000)[1]
P2000 = get_temperature_for_year(file_path, 2000)[2]

T2001 = get_temperature_for_year(file_path, 2001)[0]
H2001 = get_temperature_for_year(file_path, 2001)[1]
P2001 = get_temperature_for_year(file_path, 2001)[2]

T2002 = get_temperature_for_year(file_path, 2002)[0]
H2002 = get_temperature_for_year(file_path, 2002)[1]
P2002 = get_temperature_for_year(file_path, 2002)[2]

T2003 = get_temperature_for_year(file_path, 2003)[0]
H2003 = get_temperature_for_year(file_path, 2003)[1]
P2003 = get_temperature_for_year(file_path, 2003)[2]

T2004 = get_temperature_for_year(file_path, 2004)[0]
H2004 = get_temperature_for_year(file_path, 2004)[1]
P2004 = get_temperature_for_year(file_path, 2004)[2]

T2005 = get_temperature_for_year(file_path, 2005)[0]
H2005 = get_temperature_for_year(file_path, 2005)[1]
P2005 = get_temperature_for_year(file_path, 2005)[2]

T2006 = get_temperature_for_year(file_path, 2006)[0]
H2006 = get_temperature_for_year(file_path, 2006)[1]
P2006 = get_temperature_for_year(file_path, 2006)[2]

T2007 = get_temperature_for_year(file_path, 2007)[0]
H2007 = get_temperature_for_year(file_path, 2007)[1]
P2007 = get_temperature_for_year(file_path, 2007)[2]

T2008 = get_temperature_for_year(file_path, 2008)[0]
H2008 = get_temperature_for_year(file_path, 2008)[1]
P2008 = get_temperature_for_year(file_path, 2008)[2]

T2009 = get_temperature_for_year(file_path, 2009)[0]
H2009 = get_temperature_for_year(file_path, 2009)[1]
P2009 = get_temperature_for_year(file_path, 2009)[2]

T2010 = get_temperature_for_year(file_path, 2010)[0]
H2010 = get_temperature_for_year(file_path, 2010)[1]
P2010 = get_temperature_for_year(file_path, 2010)[2]

T2011 = get_temperature_for_year(file_path, 2011)[0]
H2011 = get_temperature_for_year(file_path, 2011)[1]
P2011 = get_temperature_for_year(file_path, 2011)[2]

T2012 = get_temperature_for_year(file_path, 2012)[0]
H2012 = get_temperature_for_year(file_path, 2012)[1]
P2012 = get_temperature_for_year(file_path, 2012)[2]

T2013 = get_temperature_for_year(file_path, 2013)[0]
H2013 = get_temperature_for_year(file_path, 2013)[1]
P2013 = get_temperature_for_year(file_path, 2013)[2]

T2014 = get_temperature_for_year(file_path, 2014)[0]
H2014 = get_temperature_for_year(file_path, 2014)[1]
P2014 = get_temperature_for_year(file_path, 2014)[2]

T2015 = get_temperature_for_year(file_path, 2015)[0]
H2015 = get_temperature_for_year(file_path, 2015)[1]
P2015 = get_temperature_for_year(file_path, 2015)[2]

T2016 = get_temperature_for_year(file_path, 2016)[0]
H2016 = get_temperature_for_year(file_path, 2016)[1]
P2016 = get_temperature_for_year(file_path, 2016)[2]

T2017 = get_temperature_for_year(file_path, 2017)[0]
H2017 = get_temperature_for_year(file_path, 2017)[1]
P2017 = get_temperature_for_year(file_path, 2017)[2]

T2018 = get_temperature_for_year(file_path, 2018)[0]
H2018 = get_temperature_for_year(file_path, 2018)[1]
P2018 = get_temperature_for_year(file_path, 2018)[2]

T2019 = get_temperature_for_year(file_path, 2019)[0]
H2019 = get_temperature_for_year(file_path, 2019)[1]
P2019 = get_temperature_for_year(file_path, 2019)[2]

T2020 = get_temperature_for_year(file_path, 2020)[0]
H2020 = get_temperature_for_year(file_path, 2020)[1]
P2020 = get_temperature_for_year(file_path, 2020)[2]

T2021 = get_temperature_for_year(file_path, 2021)[0]
H2021 = get_temperature_for_year(file_path, 2021)[1]
P2021 = get_temperature_for_year(file_path, 2021)[2]

T2022 = get_temperature_for_year(file_path, 2022)[0]
H2022 = get_temperature_for_year(file_path, 2022)[1]
P2022 = get_temperature_for_year(file_path, 2022)[2]

T2023 = get_temperature_for_year(file_path, 2023)[0]
H2023 = get_temperature_for_year(file_path, 2023)[1]
P2023 = get_temperature_for_year(file_path, 2023)[2]
print(type(P2023))

tt=0
Carying_2006=Cf2006[tt:365]
Carying_2006 = Cf2006[tt:365]
Carying_2007 = Cf2007[tt:365]
Carying_2008 = Cf2008[tt:365]
Carying_2009 = Cf2009[tt:365]
Carying_2010 = Cf2010[tt:365]
Carying_2011 = Cf2011[tt:365]
Carying_2012 = Cf2012[tt:365]
Carying_2013 = Cf2013[tt:365]
Carying_2014 = Cf2014[tt:365]
Carying_2015 = Cf2015[tt:365]
Carying_2016 = Cf2016[tt:365]
Carying_2017 = Cf2017[tt:365]
Carying_2018 = Cf2018[tt:365]
Carying_2019 = Cf2019[tt:365]
Carying_2020 = Cf2020[tt:365]

tt=0
#Carying_2000 = C2000[tt:365]
Temp_2000 = T2000[tt:365]
Humidity_2000 = H2000[tt:365]
Precipitation_2000 = P2000[tt:365]

#Carying_2001 = C2001[tt:365]
Temp_2001 = T2001[tt:365]
Humidity_2001 = H2001[tt:365]
Precipitation_2001 = P2001[tt:365]

#Carying_2002 = C2002[tt:365]
Temp_2002 = T2002[tt:365]
Humidity_2002 = H2002[tt:365]
Precipitation_2002 = P2002[tt:365]

#Carying_2003 = C2003[tt:365]
Temp_2003 = T2003[tt:365]
Humidity_2003 = H2003[tt:365]
Precipitation_2003 = P2003[tt:365]

#Carying_2004 = C2004[tt:365]
Temp_2004 = T2004[tt:365]
Humidity_2004 = H2004[tt:365]
Precipitation_2004 = P2004[tt:365]

#Carying_2005 = C2005[tt:365]
Temp_2005 = T2005[tt:365]
Humidity_2005 = H2005[tt:365]
Precipitation_2005 = P2005[tt:365]

Carying_2006=Cf2006[tt:365]
Temp_2006=T2006[tt:365]
Humidity_2006=H2006[tt:365]
Precipitation_2006=P2006[tt:365]
#----------------------------------------------------------
Carying_2007=Cf2007[tt:365]
Temp_2007=T2007[tt:365]
Humidity_2007=H2007[tt:365]
Precipitation_2007=P2007[tt:365]
#----------------------------------------------------------
Carying_2008=Cf2008[tt:365]
Temp_2008=T2008[tt:365]
Humidity_2008=H2008[tt:365]
Precipitation_2008=P2008[tt:365]
#----------------------------------------------------------
Carying_2009=Cf2009[tt:365]
Temp_2009=T2009[tt:365]
Humidity_2009=H2009[tt:365]
Precipitation_2009=P2009[tt:365]
#----------------------------------------------------------
# Carying_2010=C2010[tt:365]
# Temp_2010=T2010[tt:365]
# Humidity_2010=H2010[tt:365]
#Precipitation_2010=P2010[tt:365]
#----------------------------------------------------------
Carying_2011=Cf2011[tt:365]
Temp_2011=T2011[tt:365]
Humidity_2011=H2011[tt:365]
Precipitation_2011=P2011[tt:365]
#----------------------------------------------------------
Carying_2012=Cf2012[tt:365]
Temp_2012=T2012[tt:365]
Humidity_2012=H2012[tt:365]
Precipitation_2012=P2012[tt:365]
#----------------------------------------------------------
Carying_2013=Cf2013[tt:365]
Temp_2013=T2013[tt:365]
Humidity_2013=H2013[tt:365]
Precipitation_2013=P2013[tt:365]
#----------------------------------------------------------
Carying_2014=Cf2014[tt:365]
Temp_2014=T2014[tt:365]
Humidity_2014=H2014[tt:365]
Precipitation_2014=P2014[tt:365]
#----------------------------------------------------------
Carying_2015=Cf2015[tt:365]
Temp_2015=T2015[tt:365]
Humidity_2015=H2015[tt:365]
Precipitation_2015=P2015[tt:365]
#----------------------------------------------------------
Carying_2016=Cf2016[tt:365]
Temp_2016=T2016[tt:365]
Humidity_2016=H2016[tt:365]
Precipitation_2016=P2016[tt:365]
#----------------------------------------------------------
Carying_2017=Cf2017[tt:365]
Temp_2017=T2017[tt:365]
Humidity_2017=H2017[tt:365]
Precipitation_2017=P2017[tt:365]
#----------------------------------------------------------
Carying_2018=Cf2018[tt:365]
Temp_2018=T2018[tt:365]
Humidity_2018=H2018[tt:365]
Precipitation_2018=P2018[tt:365]
#----------------------------------------------------------
Carying_2019=Cf2019[tt:365]
Temp_2019=T2019[tt:365]
Humidity_2019=H2019[tt:365]
Precipitation_2019=P2019[tt:365]
#----------------------------------------------------------
Carying_2020=Cf2020[tt:365]
Temp_2020=T2020[tt:365]
Humidity_2020=H2020[tt:365]
Precipitation_2020=P2020[tt:365]
#----------------------------------------------------------
# Carying_2021=C2021[tt:365]
Temp_2021=T2021[tt:365]
Humidity_2021=H2021[tt:365]
Precipitation_2021=P2021[tt:365]
# #----------------------------------------------------------
# Carying_2022=C2022[tt:365]
Temp_2022=T2022[tt:365]
Humidity_2022=H2022[tt:365]
Precipitation_2022=P2022[tt:365]
#---------------------------------------------------------
Temp_2023=T2023[tt:365]
Humidity_2023=H2023[tt:365]
Precipitation_2023=P2023[tt:365]

def shifter(a,days):
    d=days
    T=a
    l=[]
    for i in range(0,365):
        if i < (365-d):
            l.append(T[i+d])
        else:
            l.append(T[i-(365-d)])
    return(l)
years = [year for year in range(2006, 2021) if year != 2010]
maxpre = []

for year in years:
    variable_name = f'P{year}'
    max_value = max(globals()[variable_name])
    maxpre.append(max_value)

#print(maxpre)

days=-15
# Precipitation_2000 = shifter(Precipitation_2000, days)
# Precipitation_2001 = shifter(Precipitation_2001, days)
# Precipitation_2002 = shifter(Precipitation_2002, days)
# Precipitation_2003 = shifter(Precipitation_2003, days)
# Precipitation_2004 = shifter(Precipitation_2004, days)
# Precipitation_2005 = shifter(Precipitation_2005, days)
Precipitation_2006=shifter(Precipitation_2006,days)
Precipitation_2007=shifter(Precipitation_2007,days)
Precipitation_2008=shifter(Precipitation_2008,days)
Precipitation_2009=shifter(Precipitation_2009,days)
#Precipitation_2010=shifter(Precipitation_2010,days)
Precipitation_2011=shifter(Precipitation_2011,days)
Precipitation_2012=shifter(Precipitation_2012,days)
Precipitation_2013=shifter(Precipitation_2013,days)
Precipitation_2014=shifter(Precipitation_2014,days)
Precipitation_2015=shifter(Precipitation_2015,days)
Precipitation_2016=shifter(Precipitation_2016,days)
Precipitation_2017=shifter(Precipitation_2017,days)
Precipitation_2018=shifter(Precipitation_2018,days)
Precipitation_2019=shifter(Precipitation_2019,days)
Precipitation_2020=shifter(Precipitation_2020,days)
Precipitation_2021=shifter(Precipitation_2021,days)
# Precipitation_2022=shifter(Precipitation_2022,days)
# Precipitation_2023=shifter(Precipitation_2023,days)

#----------------------------------------------------------
data_dict_Carying = {}
data_dict_Temp={}
data_dict_Humid={}
prerange=[i/100 for i in range(0,8000)]
for h in prerange:
    # print("h",h)
    diclist=[]
    diclisttem=[]
    diclisthum=[]
    for pre in Precipitation_2006:
        if h==pre: 
            diclist.append(Carying_2006[(Precipitation_2006.index(pre))])
            diclisttem.append(Temp_2006[(Precipitation_2006.index(pre))])
            diclisthum.append(Humidity_2006[(Precipitation_2006.index(pre))])
    for pre in Precipitation_2007:
        if h==pre: 
            diclist.append(Carying_2007[(Precipitation_2007.index(pre))])
            diclisttem.append(Temp_2007[(Precipitation_2007.index(pre))])
            diclisthum.append(Humidity_2007[(Precipitation_2007.index(pre))])
    for pre in Precipitation_2008:
        if h==pre: 
            diclist.append(Carying_2008[(Precipitation_2008.index(pre))])
            diclisttem.append(Temp_2008[(Precipitation_2008.index(pre))])
            diclisthum.append(Humidity_2008[(Precipitation_2008.index(pre))])
    for pre in Precipitation_2009:
        if h==pre: 
            diclist.append(Carying_2009[(Precipitation_2009.index(pre))])
            diclisttem.append(Temp_2009[(Precipitation_2009.index(pre))])
            diclisthum.append(Humidity_2009[(Precipitation_2009.index(pre))])
    for pre in Precipitation_2011:
        if h==pre: 
            diclist.append(Carying_2011[(Precipitation_2011.index(pre))])
            diclisttem.append(Temp_2011[(Precipitation_2011.index(pre))])
            diclisthum.append(Humidity_2011[(Precipitation_2011.index(pre))])
    for pre in Precipitation_2012:
        if h==pre: 
            diclist.append(Carying_2012[(Precipitation_2012.index(pre))])
            diclisttem.append(Temp_2012[(Precipitation_2012.index(pre))])
            diclisthum.append(Humidity_2012[(Precipitation_2012.index(pre))])
    for pre in Precipitation_2013:
        if h==pre: 
            diclist.append(Carying_2013[(Precipitation_2013.index(pre))])
            diclisttem.append(Temp_2013[(Precipitation_2013.index(pre))])
            diclisthum.append(Humidity_2013[(Precipitation_2013.index(pre))])
    for pre in Precipitation_2014:
        if h==pre: 
            diclist.append(Carying_2014[(Precipitation_2014.index(pre))])
            diclisttem.append(Temp_2014[(Precipitation_2014.index(pre))])
            diclisthum.append(Humidity_2014[(Precipitation_2014.index(pre))])
    for pre in Precipitation_2015:
        if h==pre: 
            diclist.append(Carying_2015[(Precipitation_2015.index(pre))])
            diclisttem.append(Temp_2015[(Precipitation_2015.index(pre))])
            diclisthum.append(Humidity_2015[(Precipitation_2015.index(pre))])
    for pre in Precipitation_2016:
        if h==pre: 
            diclist.append(Carying_2016[(Precipitation_2016.index(pre))])
            diclisttem.append(Temp_2016[(Precipitation_2016.index(pre))])
            diclisthum.append(Humidity_2016[(Precipitation_2016.index(pre))])
    for pre in Precipitation_2017:
        if h==pre: 
            diclist.append(Carying_2017[(Precipitation_2017.index(pre))])
            diclisttem.append(Temp_2017[(Precipitation_2017.index(pre))])
            diclisthum.append(Humidity_2017[(Precipitation_2017.index(pre))])
    for pre in Precipitation_2018:
        if h==pre: 
            diclist.append(Carying_2018[(Precipitation_2018.index(pre))])
            diclisttem.append(Temp_2018[(Precipitation_2018.index(pre))])
            diclisthum.append(Humidity_2018[(Precipitation_2018.index(pre))])
    for pre in Precipitation_2019:
        if h==pre: 
            diclist.append(Carying_2019[(Precipitation_2019.index(pre))])
            diclisttem.append(Temp_2019[(Precipitation_2019.index(pre))])
            diclisthum.append(Humidity_2019[(Precipitation_2019.index(pre))])
    for pre in Precipitation_2020:
        if h==pre: 
            diclist.append(Carying_2020[(Precipitation_2020.index(pre))])
            diclisttem.append(Temp_2020[(Precipitation_2020.index(pre))])
            diclisthum.append(Humidity_2020[(Precipitation_2020.index(pre))])
    if len(diclist)>0:
        #print(h)
        # print(diclist)
        # print(diclisttem)
        # print(diclisthum)
        data_dict_Carying[h]=list((diclist))
        data_dict_Temp[h]=list(diclisttem)
        data_dict_Humid[h]=list(diclisthum)



#---------------------------------------------------------------------------
def K_T_H(humidity,temperature,precipitation):   
    rr=precipitation
    # print("rr",rr)
    tem=temperature
    hum=humidity
    x = data_dict_Temp[rr]
    y = data_dict_Humid[rr]
    z = data_dict_Carying[rr]
    if len(x)<3:
        C=np.mean(data_dict_Carying[rr])
        # print(x)
        # print(y)
    elif len(x)>=3 : 
        # print(x)
        # print(y)
        # def polynomial_model(X, a, b, c,d,e,f,r):
        #     x, y = X
        #     return a*x**3 + b*x**2*y**2 + c*y**3 + d*x*y**2 + e*y*x**2 + f*x*y+r

        def polynomial_model(X, a, b,c):
            x, y = X

            return a*x + b*y+c
        params, params_covariance = curve_fit(polynomial_model, (x, y), z)
        x_values = np.linspace(min(x), max(x), 100)
        y_values = np.linspace(min(y), max(y), 100)
        x_mesh, y_mesh = np.meshgrid(x_values, y_values) 
        z_mesh = polynomial_model((x_mesh, y_mesh), *params)
        angles = [(30, -50), (30, 50), (90, 0)]

        # Create a 3D plot from three different angles
        # fig = plt.figure(figsize=(15, 4))

        # for i, angle in enumerate(angles, start=1):
        #     ax = fig.add_subplot(1, 3, i, projection='3d')
        #     ax.plot_surface(x_mesh, y_mesh, z_mesh, alpha=0.5, rstride=10, cstride=10, color='r', edgecolor='none')
        #     ax.scatter(x, y, z, color='b', marker='o')
        #     ax.view_init(elev=angle[0], azim=angle[1])
        #     ax.set_xlabel('Temperature')
        #     ax.set_ylabel('Humidity(g/kg)')
        #     ax.set_zlabel('Carrying capacity')
        #     ax.set_title(f'Precipitation value= {rr} mm/day')
        #     plt.savefig('predictionof carryin-orange 0.05'+'.pdf')
        # plt.show()
        C = polynomial_model((tem, hum), *params)
    #print(C)
    return(C)
#K_T_H(2.44 ,14.53, 0.05 )

#-----------------------------------function of carying capacity---------------
Year=1989
# Hum=H2023[0:365]
# dama=T2023[0:365]
# Pre=P2023[0:365]
#-------------
d0=get_temperature_for_year(file_path, Year)[0]
dama=[d0[i] for i in range(0,365)]
d1=get_temperature_for_year(file_path, Year)[1]
Hum=[d1[i] for i in range(0,365)]
d2=get_temperature_for_year(file_path, Year)[2]
Pre=[d2[i]for i in range(0,365)]

Coematrix=np.zeros([9,365]) 
    
for t in range(0, 365):
    # print(t)
    if Pre[t] > max(data_dict_Carying.keys()):
        Coematrix[8, t] = K_T_H(Hum[t], dama[t], max(data_dict_Carying.keys()))
    elif Pre[t] < min(data_dict_Carying.keys()):
        Coematrix[8, t] = K_T_H(Hum[t], dama[t], min(data_dict_Carying.keys()))
    elif Pre[t] not in data_dict_Carying.keys():
        number = Pre[t]
        my_vector = np.array(list(data_dict_Carying.keys()))  # Convert list to numpy array
        differences = np.abs(my_vector - number)
        closest_index = np.argmin(differences)
        closest_element = my_vector[closest_index]
        if K_T_H(Hum[t], dama[t], closest_element) > 0:
            Coematrix[8, t] = K_T_H(Hum[t], dama[t], closest_element)
        else:
            Coematrix[8, t] = min(data_dict_Carying[closest_element])
    else:
        if K_T_H(Hum[t], dama[t], Pre[t]) > 0:
            Coematrix[8, t] = K_T_H(Hum[t], dama[t], Pre[t])
        else:
            Coematrix[8, t] = min(data_dict_Carying[Pre[t]])
    print(t, Hum[t], dama[t], Pre[t], Coematrix[8, t])


        
def moving_average_numpy_full(data, window_size):
    """Compute the moving average using NumPy with padding."""
    window = np.ones(int(window_size)) / float(window_size)
    extended_data = np.pad(data, pad_width=(window_size//2, window_size-1-window_size//2), mode='edge')
    return np.convolve(extended_data, window, 'valid')
window_size=5
plt.plot(moving_average_numpy_full(Coematrix[8,:], window_size),'g-')
plt.plot(Coematrix[8,:],'r-')
plt.show()
Coematrix[8,:]=moving_average_numpy_full(Coematrix[8,:], window_size)
def tim2tem(time):
    tempc=dama[time]
    if tempc>34:tempc=34
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
plt.plot(R,':k',label='estimated profile of mosquitoes')

C10=Coematrix[8,:]
x=odeint(BMPPT,x0,t)
R=x[0:365,3]
I=x[0:365,2]
M=x[0:365,6]
plt.plot(I,'--k',label='estimated profile of mosquitoes')

moving_avg_values = moving_average_numpy_full(Coematrix[8, :], window_size)
moving_avg_values_column = moving_avg_values.reshape(-1, 1)
# np.savetxt(f'C_esti_{Year}.txt', moving_avg_values_column)
# np.savetxt(f'C{Year}.txt', moving_avg_values_column)






















