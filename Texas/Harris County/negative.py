import numpy as np
import pandas as pd
import statsmodels.api as sm
import matplotlib.pyplot as plt

# --------------------------
# 1. Load Historical Data (2014-2023, excluding 2020)
# --------------------------
years = [yr for yr in range(2014, 2024) if yr != 2020]  # skip 2020
hist_cases = []
hist_dates = []

for year in years:
    filename = f"./H{year}.txt"
    try:
        data = np.loadtxt(filename)
        hist_cases.extend(data)
        # Generate weekly dates for that year
        hist_dates.extend(pd.date_range(start=f"{year}-01-01", periods=len(data), freq='W'))
    except Exception as e:
        print(f"Error loading {filename}: {e}")

if len(hist_cases) == 0:
    raise ValueError("No historical data found for 2014-2023 (excluding 2020).")

hist_data = pd.DataFrame({
    'date': pd.to_datetime(hist_dates),
    'cases': np.array(hist_cases, dtype=float)
})
# Add week number and seasonal features
hist_data['week_of_year'] = hist_data['date'].dt.isocalendar().week.astype(float)
hist_data['sin_week'] = np.sin(2 * np.pi * hist_data['week_of_year'] / 52)
hist_data['cos_week'] = np.cos(2 * np.pi * hist_data['week_of_year'] / 52)
# Create a sequential time index
hist_data['time_index'] = np.arange(len(hist_data)).astype(float)

# Define a fixed scale factor from historical data.
# For example, use the maximum time_index from historical data.
scale_factor = hist_data['time_index'].max()
if scale_factor == 0:
    scale_factor = 1.0  # fallback

# Compute a scaled time variable: value between 0 and 1
hist_data['time_scaled'] = hist_data['time_index'] / scale_factor

# --------------------------
# 2. Load 2024 Data
# --------------------------
year_2024 = 2024
try:
    data_2024 = np.loadtxt(f"./H{year_2024}.txt")
except Exception as e:
    raise ValueError(f"Error loading H{year_2024}.txt: {e}")

n_weeks_2024 = len(data_2024)
dates_2024 = pd.date_range(start=f"{year_2024}-01-01", periods=n_weeks_2024, freq='W')
data_2024 = pd.DataFrame({
    'date': dates_2024,
    'cases': np.array(data_2024, dtype=float)
})
data_2024['week_of_year'] = data_2024['date'].dt.isocalendar().week.astype(float)
data_2024['sin_week'] = np.sin(2 * np.pi * data_2024['week_of_year'] / 52)
data_2024['cos_week'] = np.cos(2 * np.pi * data_2024['week_of_year'] / 52)
data_2024['time_index'] = np.arange(len(data_2024)).astype(float)
# Use the same scaling factor as for historical data:
data_2024['time_scaled'] = data_2024['time_index'] / scale_factor

# --------------------------
# 3. Rolling One-Step-Ahead Forecast for 2024
# --------------------------
# For each week t (from t=0 to n_weeks_2024-2), we use:
# training = historical data + 2024 data for weeks [0, t]
# and then predict week t+1.
predictions = []
actuals = []
forecast_dates = []

for i in range(n_weeks_2024 - 1):
    # Build training data: historical data + observed 2024 data up to week t (i)
    if i == 0:
        train_data = hist_data.copy()
    else:
        train_data = pd.concat([hist_data, data_2024.iloc[:i]], ignore_index=True)
    
    train_data = train_data.sort_values("date").reset_index(drop=True)
    
    # We'll use the fixed "time_scaled" computed with our scale_factor.
    # (For training data from historical years, it was already computed.)
    # For any new data that might not have it (should not happen), we do:
    if 'time_scaled' not in train_data.columns:
        train_data['time_scaled'] = train_data['time_index'] / scale_factor

    # Build design matrix and response.
    # Use predictors: intercept, time_scaled, sin_week, cos_week.
    X_train = train_data[['time_scaled', 'sin_week', 'cos_week']]
    X_train = sm.add_constant(X_train)
    y_train = train_data['cases']
    
    # Ensure numeric types
    X_train = X_train.astype(float)
    y_train = y_train.astype(float)
    
    # Fit Negative Binomial model using GLM
    nb_model = sm.GLM(y_train, X_train, family=sm.families.NegativeBinomial()).fit()
    
    # Forecast for week t+1: use row at index i+1 from 2024 data
    forecast_row = data_2024.iloc[i+1]
    # The forecast time_scaled is the value in the 2024 data
    ts_new = forecast_row['time_scaled']
    X_new = pd.DataFrame({
        'const': [1.0],
        'time_scaled': [ts_new],
        'sin_week': [np.sin(2 * np.pi * forecast_row['week_of_year'] / 52)],
        'cos_week': [np.cos(2 * np.pi * forecast_row['week_of_year'] / 52)]
    })
    X_new = X_new.astype(float)
    
    # Predict the mean count using the log link (exp(linear predictor))
    mu_pred = nb_model.predict(X_new)[0]
    
    # Optionally round to nearest integer
    predictions.append(int(round(mu_pred)))
    actuals.append(int(round(forecast_row['cases'])))
    forecast_dates.append(forecast_row['date'])

# --------------------------
# 4. Plot Forecasts vs. Actual Values
# --------------------------
plt.figure(figsize=(10, 6))
plt.plot(forecast_dates, predictions, marker='o', label='Predicted Cases')
plt.plot(forecast_dates, actuals, marker='x', linestyle='--', label='Actual Cases')
plt.xlabel('Date')
plt.ylabel('Weekly Cases')
plt.title('Rolling One-Step-Ahead Negative Binomial Forecasts for 2024')
plt.legend()
plt.grid(True)
plt.xticks(rotation=45)
plt.tight_layout()
plt.show()

print("Predicted weekly cases for 2024 (as integers):")
print(len(predictions))
