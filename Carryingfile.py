




import numpy as np
from scipy.integrate import odeint
from scipy.optimize import minimize_scalar
import matplotlib.pyplot as plt
import pandas as pd
import os
from scipy.interpolate import interp1d




def moving_average(a, window=7):
    kernel = np.ones(window) / window
    return np.convolve(a, kernel, mode='same')

YEARS = [
    2006,2007,2008,2009,2010,2011,2012,2013,
    2015,2016,2017,2018,2019,2020,2021,2022,2023,2024
]

control=0
controlday=0
all_C_ens = []  # Store ensemble carrying capacity for each year

climate_df = pd.read_csv('orange.csv')

def temperature_series(year):
    t = climate_df.loc[climate_df['Year'] == year, 'T'].to_numpy()
    if t.size == 0:
        raise ValueError(f'Year {year} not found in CSV')
    return t

def build_coef(T_series):
    def fwrap(fun):
        return np.vectorize(lambda d: fun(T_series[int(d)]))
    gamma  = fwrap(lambda T: 0 if (T<=4.3 or T>=39.9) else T*4.12e-5*(T-4.3)*np.sqrt(39.9-T))
    pLST   = fwrap(lambda T: 0 if (T<=5.9 or T>=43.1) else -2.12e-3*(T-5.9)*(T-43.1))
    muM    = fwrap(lambda T: 1/1e-16 if T>41.3 else (1/45 if T<14 else 1/(0.99*(-1.69*T+69.6))))
    BRT    = fwrap(lambda T: 0 if (T<=2.3 or T>=32) else T*1.67e-4*(T-2.3)*np.sqrt(32-T))
    alpha  = fwrap(lambda T: 0 if (T<=5.3 or T>=38.9) else -5.98e-1*(T-5.3)*(T-38.9)*BRT(int(T)))
    bET    = fwrap(lambda T: 0 if (T<=3.2 or T>=42.6) else -2.11e-3*(T-3.2)*(T-42.6))
    PDRT   = fwrap(lambda T: 0 if (T<=11.2 or T>=44.7) else T*6.57e-5*(T-11.2)*np.sqrt(44.7-T))
    VCT    = fwrap(lambda T: 0 if (T<=11.3 or T>=41.9) else -2.94e-3*(T-11.3)*(T-41.9))
    HR     = lambda d: 1.0
    muegg  = lambda d: 0.0
    ADR    = lambda d: pLST(d)*gamma(d)
    muAQ   = lambda d: gamma(d)*(1-pLST(d))

    funcs = (gamma, pLST, muM, BRT, alpha, bET,
             PDRT, VCT, HR, muegg, ADR, muAQ)
    M = np.zeros((12, 365))
    for day in range(365):
        M[:, day] = [f(day) if not callable(f) else f(day) for f in funcs]
    return M

def rhs_factory(C_ac_vec, controlday=0):
    def BMPPT(x, t):
        day = int(t) % 365
        C_ac = C_ac_vec[int(t)]
        muB = 0.15/365; muhum = 0.01; incH = 1/6; infH = 0.14; muH = 1/(77*365)
        incB = 1/3; recC = 1/6; muwnvB = 0.9
        SH, EH, IH, RH, EM, AM, MS, ME, MI, BS, BE, BI, BR = x
        TH = SH+EH+IH+RH; D = BS+BE+BI+BR; eps = 1e-8
        ph = (Coematrix[3,day]*Coematrix[7,day])/(D+TH+eps)
        dsH = muH*TH - ph*MI*SH - muH*SH
        deH = 0 if day<controlday else ph*MI*SH - incH*EH - muH*EH
        diH = 0 if day<control else incH*EH - infH*IH - muhum*IH - muH*IH
        drH = infH*IH - muH*RH
        degg = Coematrix[4,day]*(MS+ME+MI) - Coematrix[5,day]*EM \
               - Coematrix[8,day]*EM - Coematrix[9,day]*EM
        dacu = Coematrix[5,day]*EM*max(0,1-AM/(C_ac+eps)) \
               - Coematrix[0,day]*Coematrix[1,day]*AM \
               - Coematrix[10,day]*AM - Coematrix[11,day]*AM
        dsM  = Coematrix[10,day]*AM - ph*BI*MS - Coematrix[2,day]*MS
        deM  = ph*BI*MS - Coematrix[6,day]*ME - Coematrix[2,day]*ME
        diM  = Coematrix[6,day]*ME - Coematrix[2,day]*MI
        dsB  = muB*(BS+BE+BI+BR) - ph*Coematrix[7,day]*MI*BS
        deB  = ph*Coematrix[7,day]*MI*BS - incB*BE - muB*BE
        diB  = incB*BE - muwnvB*BI - muB*BI - recC*BI
        drB  = recC*BI - muB*BR
        return np.array([dsH,deH,diH,drH,degg,dacu,dsM,deM,diM,dsB,deB,diB,drB])
    return BMPPT

def mse_for_k(k0, ens_size=300, seed=0):
    rng = np.random.default_rng(seed)
    n_days = 364
    x0 = np.array([3e5,0,0,0,1,1,100,1,1,1e4,1,1,1])
    C_ac_init = np.abs(k0 + 1e1*rng.standard_normal(ens_size))
    states    = np.zeros((ens_size, n_days, 13))
    C_ens     = np.zeros((ens_size, n_days))
    states[:,0,:] = x0
    C_ens[:,0]   = C_ac_init
    obs_R = 2.0
    for t in range(1, n_days):
        for n in range(ens_size):
            C_ens[n,t] = np.abs(C_ens[n,t-1] + 5e5*rng.standard_normal())
            rhs = rhs_factory(np.full(n_days, C_ens[n,t]))
            states[n,t,:] = np.maximum(states[n,t-1,:] + rhs(states[n,t-1,:], t-1)
                                        + 1e-3*rng.standard_normal(13), 0)
        RH_pred = states[:,t,3]
        y = recovered_obs[t] if t < len(recovered_obs) else recovered_obs[-1]
        X = np.hstack([states[:,t,:], C_ens[:,t:t+1]])
        HX = RH_pred.reshape(-1,1)
        Xmean = X.mean(axis=0); HXmean = HX.mean(axis=0)
        PfHT = ((X - Xmean).T @ (HX - HXmean)) / (ens_size-1)
        HPfH = HX.var(ddof=1) + obs_R
        K = PfHT / HPfH
        innov = y + np.sqrt(obs_R)*rng.standard_normal(ens_size) - RH_pred
        X += (K @ innov.reshape(1,-1)).T
        states[:,t,:] = np.maximum(X[:,:-1], 0)
        C_ens[:,t]    = np.abs(X[:,-1])
    # deterministic ODE with mean C_ac
    C_mean = C_ens.mean(axis=0)
    rhs_det = rhs_factory(C_mean)
    x_det = odeint(rhs_det, x0, np.arange(n_days))
    RH_det = x_det[:,3]
    return np.sqrt(np.sum((RH_det - recovered_obs)**2)), C_ens

# AA = []
# BB = []

# for Year in YEARS:
#     dama = temperature_series(Year)
#     weekly_cases = np.loadtxt(f'./orangecases{Year}.txt')
#     daily_inc = np.zeros(364)
#     for wk, val in enumerate(weekly_cases):
#         if wk*7 < len(daily_inc):
#             daily_inc[wk*7] = val
#     cases = daily_inc
#     idx = [i for i in range(0,52)]
#     recovered_obs = np.cumsum(cases)
#     Coematrix = build_coef(dama)

#     # Fast 1-D optimisation on log10(k0)
#     log_lo, log_hi = 6.0, 8.5
#     obj = lambda logk: mse_for_k(10**logk, ens_size=300)[0]

#     result = minimize_scalar(obj, bounds=(log_lo, log_hi), method='bounded',
#                               options=dict(maxiter=30))
#     best_k0 = 10**result.x
#     print(f"Fast search ⇒ best k0 ≈ {best_k0:,.0f}  (MSE ≈ {result.fun:.2f})")

#     # Final run with full 100-member EnKF and get ensemble
#     final_mse, C_ens = mse_for_k(best_k0, ens_size=300, seed=123)
#     print(f"Final MSE with 100-member EnKF = {final_mse:.2f}")

#     rhs_final = rhs_factory(np.full(364, best_k0))
#     x_det = odeint(rhs_final, np.array([3e5,0,0,0,1,1,100,1,1,1e4,1,1,1]), np.arange(364))

#     plt.figure(figsize=(8,4))
#     plt.plot(recovered_obs, 'ro', ms=3, label='Observed cumulative RH')
#     plt.plot(x_det[:,3], 'k-', lw=2, label='ODE RH (best $k_0$)')
#     plt.xlabel('Day'); plt.ylabel('Cumulative recovered')
#     plt.title(f'Best initial carrying capacity = {best_k0:,.0f} - {Year}')
#     plt.legend()
#     plt.tight_layout()
#     plt.show()

#     case_days = np.where(cases > 0)[0]
#     max_case_day = np.argmax(cases)
#     for d in case_days:
#         col = 'red' if d == max_case_day else 'k'
#         alp = 0.8 if d == max_case_day else 0.3
#         lw = 1.2 if d == max_case_day else 0.8
#         plt.axvline(d, color=col, alpha=alp, lw=lw)
#     plt.plot(x_det[:,6]+x_det[:,7]+x_det[:,8], 'b-', lw=2, label=f'MT (best $k_0$) - {Year}')
#     plt.xlabel('Day'); plt.ylabel('MT')
#     plt.legend()
#     plt.tight_layout()
#     plt.show()

#     #week_idx = np.arange(30, 51)
#     cum_M = x_det[:,6]+x_det[:,7]+x_det[:,8]
#     cum_M = moving_average(cum_M, window=14)

#     # AA.append(cum_M[week_idx * 7])
#     # BB.append(weekly_cases[week_idx])

#     # print("year", Year)
#     # print("weeks with nonzero cases", idx)
#     # print("Cum M", len(cum_M[week_idx * 7]))
#     # print("Cases", len(weekly_cases[week_idx]), '\n')

#     all_C_ens.append(C_ens)  # Save ensemble for this year

#     # --- Save ensemble carrying capacity with temperature ---
#     temp_rounded = np.round(dama[:364], 1)  # Round temperature to 1 decimal place
#     save_array = np.column_stack([temp_rounded, C_ens.T])  # Shape: (n_days, n_ensemble+1)
#     np.savetxt(f'C_ac_ensemble_{Year}.txt', save_array, delimiter=',',
#                 header=','.join(['Temp'] + [f'Ens_{i}' for i in range(C_ens.shape[0])]),
#                 comments='', fmt='%.1f' + ',%.6e' * C_ens.shape[0])

# # --- Plot ensemble of carrying capacity trajectories for all years ---
# plt.figure(figsize=(14, 8))
# for i, year in enumerate(YEARS):
#     if i >= len(all_C_ens):
#         break
#     C_ens = all_C_ens[i]  # shape: (ensemble_size, n_days)
#     n_ens, n_days = C_ens.shape
#     days = np.arange(n_days)
#     # Plot all ensemble members for year with transparency
#     for ens_idx in range(n_ens):
#         plt.plot(days, C_ens[ens_idx], color='blue', alpha=0.05)
#     # Plot ensemble mean for year
#     plt.plot(days, C_ens.mean(axis=0), label=str(year), lw=2)

# plt.xlabel('Day of year')
# plt.ylabel('Carrying capacity (C_ac)')
# plt.title('Kalman Filter Estimated Ensemble of Carrying Capacity (best $k_0$)')
# plt.legend()
# plt.tight_layout()
# plt.show()

# ----------- 1. Read and Stack All Ensembles by Temperature -----------

all_T = []
all_C = []
for Y in YEARS:
    fname = f'C_ac_ensemble_{Y}.txt'
    if not os.path.exists(fname):
        print(f"File {fname} not found, skipping.")
        continue
    arr = np.loadtxt(fname, delimiter=',', skiprows=1)  # Skip header row
    # First column: temperature; rest: ensemble carrying capacities
    T_this = arr[:, 0]        # (n_days,)
    C_this = arr[:, 1:]       # (n_days, n_ensemble)
    # Repeat temperature values to match each ensemble member per day
    T_rep = np.repeat(T_this, C_this.shape[1])
    C_flat = C_this.flatten()
    all_T.append(T_rep)
    all_C.append(C_flat)

# Concatenate all years
all_T = np.concatenate(all_T)
all_C = np.concatenate(all_C)

# ----------- 2. Bin by Temperature (0.1°C resolution) -----------

bin_edges = np.arange(np.floor(all_T.min()), np.ceil(all_T.max()) + 0.1, 0.1)
bin_centers = (bin_edges[:-1] + bin_edges[1:]) / 2

T_bin_func = []
mean_func = []
q25_func = []
median_func = []
q75_func = []

# ----------- 3. Calculate Quantiles in Each Temperature Bin -----------

for i, (t_low, t_high) in enumerate(zip(bin_edges[:-1], bin_edges[1:])):
    mask = (all_T >= t_low) & (all_T < t_high)
    C_in_bin = all_C[mask]
    if len(C_in_bin) < 5:
        continue  # Skip bins with too few samples
    T_bin_func.append(bin_centers[i])
    mean_func.append(np.mean(C_in_bin))
    q25_func.append(np.percentile(C_in_bin, 25))
    median_func.append(np.percentile(C_in_bin, 50))
    q75_func.append(np.percentile(C_in_bin, 80))

# ----------- 4. Smooth Quantile Curves -----------

def moving_average(arr, window_size=5):
    return pd.Series(arr).rolling(window=window_size, center=True, min_periods=1).mean().values

window_size = 30
mean_func = moving_average(mean_func, window_size)
q25_func = moving_average(q25_func, window_size)
median_func = moving_average(median_func, window_size)
q75_func = moving_average(q75_func, window_size)

# ----------- 5. Plot Carrying Capacity as Function of Temperature -----------

plt.figure(figsize=(10, 6))
plt.plot(T_bin_func, mean_func, 'm-', label='Mean', lw=2)
plt.plot(T_bin_func, median_func, 'r-', label='Median', lw=2)
plt.plot(T_bin_func, q25_func, 'b--', label='Q25', lw=1.5)
plt.plot(T_bin_func, q75_func, 'g--', label='Q75', lw=1.5)
plt.fill_between(T_bin_func, q25_func, q75_func, color='gray', alpha=0.25, label='IQR (Q25-Q75)')
plt.xlabel("Temperature (°C)")
plt.ylabel("Carrying Capacity $C_{ac}$")
plt.title("Carrying Capacity $C_{ac}$ as a Function of Temperature")
plt.legend()
plt.grid(True)
plt.tight_layout()
plt.show()






















# Choose which estimate to use for carrying capacity function
# Options: mean_func, median_func, q25_func, q75_func
# Example uses median_func here:

Cac_func_values =q75_func# median_func  # as you specified
T_func_values = np.array(T_bin_func)

# Stack temperature and carrying capacity into columns
save_array = np.column_stack([T_func_values, Cac_func_values])

# Save to text file with header
np.savetxt('C_predict_80.txt', save_array, delimiter=',', 
           header='Temperature,CarryingCapacity', comments='', fmt='%.2f,%.6e')




# Interpolate the carrying capacity function for any temperature
Cac_interp = interp1d(T_func_values, Cac_func_values,
                      kind='linear', bounds_error=False,
                      fill_value=(Cac_func_values[0], Cac_func_values[-1]))

def carrying_capacity_from_temp(T_daily):
    """Given daily temperature array, return corresponding daily carrying capacity C_ac."""
    return Cac_interp(T_daily)

# Run model for all years using temperature-dependent carrying capacity function
for Year in YEARS:
    print(f"Running model for Year {Year} with temp-dependent carrying capacity...")

    dama = temperature_series(Year)
    weekly_cases = np.loadtxt(f'./orangecases{Year}.txt')
    daily_inc = np.zeros(364)
    for wk, val in enumerate(weekly_cases):
        if wk*7 < len(daily_inc):
            daily_inc[wk*7] = val
    recovered_obs = np.cumsum(daily_inc)

    # Calculate carrying capacity C_ac daily from temperature
    C_ac_daily = carrying_capacity_from_temp(dama[:364])

    # Build coefficient matrix
    Coematrix = build_coef(dama)

    # RHS function with daily C_ac from temperature
    rhs = rhs_factory(C_ac_daily)

    # Initial state vector
    x0 = np.array([3e5, 0, 0, 0, 1, 1, 100, 1, 1, 1e4, 1, 1, 1])

    # Solve ODE system
    t_span = np.arange(364)
    x_sol = odeint(rhs, x0, t_span)

    RH_model = x_sol[:, 3]  # Recovered compartment

    # Plot observed vs modeled recovered cases
    plt.figure(figsize=(10, 5))
    plt.plot(recovered_obs, 'ro', label='Observed cumulative recovered')
    plt.plot(RH_model, 'b-', lw=2, label='Model cumulative recovered (temp-dependent $C_{ac}$)')
    plt.xlabel('Day of Year')
    plt.ylabel('Cumulative Recovered')
    plt.title(f'Year {Year}: Model Fit Using Temp-Dependent Carrying Capacity')
    plt.legend()
    plt.tight_layout()
    plt.show()

    # ----------------- ADDITION: Plot Mosquito Abundance -----------------
    # ----------------- ADDITION: Smooth, Save and Plot Mosquito Abundance -----------------
    Mosquito_abundance = x_sol[:, 6] + x_sol[:, 7] + x_sol[:, 8]  # MS + ME + MI

    # Smooth mosquito abundance with moving average (window=14)
    def moving_average(a, window=14):
        kernel = np.ones(window) / window
        return np.convolve(a, kernel, mode='same')
    
    Mosquito_abundance_smoothed = moving_average(Mosquito_abundance, window=1)

    # Save smoothed mosquito abundance to file
    #np.savetxt(f'MOR_predicted_{Year}.txt', Mosquito_abundance_smoothed, fmt='%.6e')

    # Plot smoothed mosquito abundance
    plt.figure(figsize=(10, 5))
    plt.plot(Mosquito_abundance_smoothed, 'g-', lw=2, label='Smoothed Mosquito Abundance (MS+ME+MI)')
    plt.xlabel('Day of Year')
    plt.ylabel('Mosquito Abundance')
    plt.title(f'Year {Year}: Smoothed Mosquito Abundance Based on Temp-Dependent Carrying Capacity')
    plt.legend()
    plt.tight_layout()
    plt.show()




















