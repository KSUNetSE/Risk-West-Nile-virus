import numpy as np
import pandas as pd
import statsmodels.api as sm
from statsmodels.genmod.families import NegativeBinomial

# --- Load Historical Data (2014-2023) ---
hist_cases = []
hist_dates = []

for year in range(2014, 2024):  # 2014 to 2023
    filename = f"./D{year}.txt"
    try:
        year_data = np.loadtxt(filename)
        hist_cases.extend(year_data)
        # Generate weekly dates (frequency 'W')
        hist_dates.extend(pd.date_range(start=f"{year}-01-01", periods=len(year_data), freq='W'))
    except Exception as e:
        print(f"Error loading {filename}: {e}")

if not hist_cases:
    raise ValueError("No historical data loaded from 2014-2023.")

hist_data = pd.DataFrame({
    'date': pd.to_datetime(hist_dates),
    'cases': np.array(hist_cases)
})
hist_data['week_of_year'] = hist_data['date'].dt.isocalendar().week
hist_data['sin_week'] = np.sin(2 * np.pi * hist_data['week_of_year'] / 52)
hist_data['cos_week'] = np.cos(2 * np.pi * hist_data['week_of_year'] / 52)

# --- Load 2024 Data (Full Year) ---
year_2024 = 2024
data_2024_array = np.loadtxt(f"./D{year_2024}.txt")
n_weeks = len(data_2024_array)  # expected to be 52 weeks
dates_2024 = pd.date_range(start=f"{year_2024}-01-01", periods=n_weeks, freq='W')

data_2024 = pd.DataFrame({
    'date': dates_2024,
    'cases': data_2024_array
})
data_2024['week_of_year'] = data_2024['date'].dt.isocalendar().week
data_2024['sin_week'] = np.sin(2 * np.pi * data_2024['week_of_year'] / 52)
data_2024['cos_week'] = np.cos(2 * np.pi * data_2024['week_of_year'] / 52)

# --- Rolling Forecast for Each Week of 2024 ---
predictions = []  # List to store predicted values (length 52)

for i in range(n_weeks):
    # Use historical data plus observed weeks of 2024 up to week i
    if i == 0:
        train_data = hist_data.copy()
    else:
        train_data = pd.concat([hist_data, data_2024.iloc[:i]], ignore_index=True)
    
    train_data = train_data.sort_values('date').reset_index(drop=True)
    # Create a sequential time index and scale it
    train_data['time_index'] = np.arange(len(train_data))
    train_data['time_index_scaled'] = (train_data['time_index'] - train_data['time_index'].mean()) / train_data['time_index'].std()
    
    # Prepare predictors and response for model fitting
    X_train = train_data[['time_index_scaled', 'sin_week', 'cos_week']]
    X_train = sm.add_constant(X_train)
    y_train = train_data['cases']
    X_train = X_train.astype(float)
    y_train = y_train.astype(float)
    
    # Fit GLM Negative Binomial model
    model = sm.GLM(y_train, X_train, family=NegativeBinomial())
    nb_model = model.fit()
    
    # Forecast for week i (the next week to be predicted)
    target = data_2024.iloc[i]
    
    # The forecast time index is the next index after training data
    forecast_index = len(train_data)
    forecast_index_scaled = (forecast_index - train_data['time_index'].mean()) / train_data['time_index'].std()
    
    # Use the week number from the target to compute seasonal features
    week = target['week_of_year']
    new_X = pd.DataFrame({
        'const': [1],
        'time_index_scaled': [forecast_index_scaled],
        'sin_week': [np.sin(2 * np.pi * week / 52)],
        'cos_week': [np.cos(2 * np.pi * week / 52)]
    })
    new_X = new_X.astype(float)
    new_X = np.asarray(new_X, dtype=float)
    
    # Make prediction and round to nearest nonnegative integer
    pred = nb_model.predict(new_X)[0]
    predictions.append(max(0, int(round(pred))))

# Print the list of predicted values for 52 weeks
print("Predicted weekly cases for 2024 (list of 52 values):")
print(predictions)
