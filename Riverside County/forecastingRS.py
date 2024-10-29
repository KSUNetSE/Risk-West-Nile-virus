# -*- coding: utf-8 -*-
"""
Created on Tue Apr 23 13:50:47 2024

@author: shosseini
"""

import numpy as np
import pandas as pd
import statsmodels.api as sm
from statsmodels.tsa.ar_model import AutoReg
import matplotlib.pyplot as plt

year=2023#the target year tp be predicted
def get_temperature_for_year(file_path, year):
    # Read the CSV file
    df = pd.read_csv(file_path)
    df_year = df[df['Year'] == year]
    # Convert the columns to lists of real values (floats)
    t = df_year['T'].tolist()
    h = df_year['H'].tolist()
    p = df_year['P'].tolist()
    return t, h, p

def get_data_from_years(file_path, year):
    temperature_list = []
    humidity_list = []
    precipitation_list = []

    for y in range(year-4, year):
        t, h, p = get_temperature_for_year(file_path, y)
        temperature_list.extend(t)
        humidity_list.extend(h)
        precipitation_list.extend(p)
    
    return temperature_list, humidity_list, precipitation_list
window_size=3
file_path = 'Riversideinfo.csv'
# Load the dataset

T=get_data_from_years(file_path, year)[0]
H=get_data_from_years(file_path, year)[1]
P=get_data_from_years(file_path, year)[2]
data=T 
#----------------------------------------------tem
# Fit an AR model with a lag of 365 days
model = AutoReg(data, lags=365)
model_fit = model.fit()

# Make forecast for the next 365 days
forecast = model_fit.predict(start=len(data), end=len(data)+364, dynamic=True)

# plt.savefig('Orabnge tem forecast2023'+'.pdf')
# plt.savefig('Orabnge tem forecast2021'+'.pdf')
#plt.savefig('Orabnge tem forecast2021'+'.pdf')
def moving_average_numpy_full(data, window_size):
    """Compute the moving average using NumPy with padding."""
    window = np.ones(int(window_size)) / float(window_size)
    extended_data = np.pad(data, pad_width=(window_size//2, window_size-1-window_size//2), mode='edge')
    return np.convolve(extended_data, window, 'valid')

forecast=moving_average_numpy_full(forecast,window_size)
plt.plot(forecast,'g')
plt.xlabel('Day')
plt.ylabel('Temperature')
plt.title('Forecast of Temperature for 2023 in RiverSide Counti, California')
plt.show()
with open('2023forecasted_Tem.txt', 'w') as file:
    # Append a newline to each string in the list
    file.writelines(f"{item}\n" for item in forecast)
#-------------------------------------------------Humidity
dataH=H 
# Fit an AR model with a lag of 365 days
modelH = AutoReg(dataH, lags=365)
model_fit = modelH.fit()

# Make forecast for the next 365 days
forecast = model_fit.predict(start=len(data), end=len(data)+364, dynamic=True)

# plt.savefig('Orabnge tem forecast2021'+'.pdf')
# plt.savefig('Orabnge tem forecast2021'+'.pdf')
#plt.savefig('Orabnge tem forecast2021'+'.pdf')
forecast=moving_average_numpy_full(forecast,window_size)
plt.plot(forecast,'g')
plt.xlabel('Day')
plt.ylabel('Humidity')
plt.title('Forecast of Humidity for 2023 in RiverSide Counti, California')
plt.show()
with open('2023forecasted_Hum.txt', 'w') as file:
    # Append a newline to each string in the list
    file.writelines(f"{item}\n" for item in forecast)
#--------------------------------------------------------precipitation
dataP=P
# Fit an AR model with a lag of 365 days
modelP = AutoReg(dataP, lags=365)
model_fit = modelP.fit()

# Make forecast for the next 365 days
forecast = model_fit.predict(start=len(data), end=len(data)+364, dynamic=True)

# plt.savefig('Orabnge tem forecast2021'+'.pdf')
# plt.savefig('Orabnge tem forecast2021'+'.pdf')
#plt.savefig('Orabnge tem forecast2021'+'.pdf')

forecast=moving_average_numpy_full(forecast,window_size)
plt.plot(forecast,'g')
plt.xlabel('Day')
plt.ylabel('Precipitation')
plt.title('Forecast of Precipitation for 2023 in RiverSide Counti, California')
plt.show()
with open('2023forecasted_Pre.txt', 'w') as file:
    # Append a newline to each string in the list
    file.writelines(f"{item}\n" for item in forecast)


