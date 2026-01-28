
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from scipy.stats import gaussian_kde
from scipy.integrate import odeint
from scipy.interpolate import interp1d

# 1. LOAD ALL DATA
climate_df = pd.read_csv('orange.csv')  # Contains columns: 'Year' and 'T' (temperature)
C_predict_data = np.loadtxt('C_predict.txt', delimiter=',', skiprows=1)  # Temp vs carrying capacity table
T_vals = C_predict_data[:, 0]  # temperature grid
C_vals = C_predict_data[:, 1]  # corresponding carrying capacity values
# Create an interpolation function so we can estimate carrying capacity for *any* temperature
Cac_interp = interp1d(T_vals, C_vals, kind='linear', bounds_error=False, fill_value='extrapolate')

# Define which years to use
YEARS = [
    2006, 2007, 2008, 2009, 2010, 2011, 2012, 2013,2014,
    2015, 2016, 2017, 2018, 2019, 2020,  2021,2022
    ,2023
    ,2024
    
]



case_first = {}          # year  →  first-case week (None if no cases)
mi=[]
for yr in YEARS:
    fname = list(np.loadtxt(f'orangecases{yr}.txt'))
    mi.append(7*(((np.nonzero(fname)[0])[0])))
#print("min index",min(mi))


min_index=min(mi)


# 2. HELPER FUNCTIONS
def contiguous_blocks(idxs):
    """
    Split a 1-D sorted integer array into sub-arrays of consecutive values.
    Example: [10,11,12,19,20] → [[10,11,12], [19,20]]
    """
    if idxs.size == 0:
        return []
    # positions where the gap > 1
    cuts = np.where(np.diff(idxs) > 1)[0] + 1
    return np.split(idxs, cuts)



def temperature_series(year):
    # Returns the daily temperature series for a specific year
    temps = climate_df.loc[climate_df['Year'] == year, 'T'].to_numpy()
    if temps.size == 0:
        raise ValueError(f"No temperature data for Year={year}")
    return temps

def first_case_week(year: int):
    # Loads weekly human case data for the given year and returns the week index of the first human case.
    wk = np.loadtxt(f'orangecases{year}.txt')
    if len(wk) < 52:
        wk = np.pad(wk, (0, 52 - len(wk)))
    nz = np.nonzero(wk)[0]
    return None if nz.size == 0 else nz[0]

# def average_temperature_excluding_year(excluded_year):
#     # Computes the average daily temperature across all years *except* the test year.
#     temps_list = []
#     for yr in YEARS:
#         if yr != excluded_year:
#             temps = temperature_series(yr)
#             temps_list.append(temps[:364])  # Ensures the same length for all years
#     avg_temp = np.mean(np.vstack(temps_list), axis=0)
#     return avg_temp



# def average_temperature_excluding_year(excluded_year, window=5):
#     """
#     Mean daily temperature over the `window` years that
#     immediately precede `excluded_year` (e.g. for 2006 → 2001-2005).

#     Raises an error only if temperature data for any of those
#     years is missing.
#     """
#     prev_years = [excluded_year - i for i in range(1, window + 1)]

#     temps_list = []
#     for yr in prev_years:
#         try:
#             temps = temperature_series(yr)[:364]   # ensure same length
#         except ValueError:
#             raise ValueError(f"No temperature data for Year={yr}, "
#                               f"needed to compute 5-year average before {excluded_year}.")
#         temps_list.append(temps)

#     return np.mean(np.vstack(temps_list), axis=0)

def average_temperature_excluding_year(excluded_year):
    window=5
    """
    Return a 364-day temperature series for `excluded_year` in which
    the first 150 days are the observed values for that year and the
    remaining days are replaced by the daily mean of the `window`
    years immediately preceding it.

    Parameters
    ----------
    excluded_year : int
        The year being cross-validated (its first 150 days are kept).
    window : int, default 5
        Number of prior years used to compute the climatological tail.

    Returns
    -------
    np.ndarray
        Length-364 array of daily temperatures.
    """
    # ----- 1. keep the first 150 days of the excluded year --------------------
    full_year = temperature_series(excluded_year)[:364]   # ensure 364-day length
    head = full_year[:140]

    # ----- 2. build the 5-year climatological mean for the remaining days ----
    prev_years = [excluded_year - i for i in range(1, window + 1)]
    temps_list = []
    for yr in prev_years:
        try:
            temps_list.append(temperature_series(yr)[:364])
        except ValueError:
            raise ValueError(f"Missing temperature data for year {yr}, "
                              f"needed for climatology before {excluded_year}.")
    clim_mean = np.mean(np.vstack(temps_list), axis=0)
    tail = clim_mean[140:]

    # ----- 3. concatenate head + tail ----------------------------------------
    return np.concatenate([head, tail])

def select_block_or_all(arr, x):
    """
    Parameters
    ----------
    arr : list[int]  (must be sorted ascending)
    x   : int        value to look for

    Returns
    -------
    numpy.ndarray
        • If `x` lies inside one of the contiguous runs of `arr`,
          that run is returned.
        • Otherwise the original array is returned unchanged.
    """
    tr=0
    arr = np.asarray(arr)
    if arr.size == 0:
        return arr, tr               # nothing to check

    # --- split arr into contiguous runs -----------------------------
    split_idx = np.where(np.diff(arr) > 1)[0] + 1
    blocks    = np.split(arr, split_idx)

    # --- find the block containing x (if any) -----------------------
    for blk in blocks:
        if blk[0] <= x <= blk[-1]:   # x is inside this run
            #print("choose",blk)
            tr=1
            #print("tr",tr)
            return blk, tr
        

    # x wasn’t in any block → return the whole list
    return (arr, tr)



# 3. MODEL FUNCTIONS

def build_coef(T_series):
    # Builds a matrix of model coefficients that depend on temperature for each day.
    def fwrap(fun):
        return np.vectorize(lambda d: fun(T_series[int(d)]))
    # Define all the biological/epidemiological rates as functions of temperature
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
    funcs = (gamma, pLST, muM, BRT, alpha, bET, PDRT, VCT, HR, muegg, ADR, muAQ)
    M = np.zeros((12, 365))
    n_days = len(T_series)
    for day in range(n_days):
        M[:, day] = [f(day) if not callable(f) else f(day) for f in funcs]
    return M

def rhs_factory(C_ac_vec, Coematrix, control=0, controlday=0):
    # Constructs the right-hand-side for the ODE system, using temperature-dependent parameters and carrying capacity.
    def BMPPT(x, t):
        day = int(t) % 365
        idx = min(int(t), len(C_ac_vec)-1)  # avoids index out of bounds
        C_ac = C_ac_vec[idx]
        # Model parameters for all compartments (humans, birds, mosquitoes)
        muB = 0.15/365; muhum = 0.01; incH = 1/6; infH = 0.14; muH = 1/(77*365)
        incB = 1/3; recC = 1/6; muwnvB = 0.9
        SH, EH, IH, RH, EM, AM, MS, ME, MI, BS, BE, BI, BR = x
        TH = SH+EH+IH+RH; D = BS+BE+BI+BR; eps = 1e-8
        ph = (Coematrix[3,day]*Coematrix[7,day])/(D+TH+eps)
        # Equations for each compartment
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

# 4. MAIN CROSS-VALIDATION LOOP

u = 0  # Counter for years where the method "works"

information=[]
# You can set a specific range of years for testing if desired.














for test_year in [2024]:#YEARS:
    #print(f"\n--- ------------------------------------------------Test Year {test_year} ---")
    M_train, T_train = [], []
    w0_test = None  # Index of first case for the test year
    for yr in YEARS:
        w0 = first_case_week(yr)
        if w0 is None:
            continue
        if yr == test_year:
            w0_test = w0
        else:
            #print("year = ",test_year,yr)
            # Collect 4 days from week before of vreporting
            day_end = w0 * 7 -7   # day of first case
            day_start = max(0, day_end-10)
            M = np.loadtxt(f'MOR_{yr}.txt')  # mosquito daily abundance
            T_daily = temperature_series(yr)
            T = T_daily[:len(M)]
            days = np.arange(day_start, day_end + 1)
            M10 = M[days]
            T10 = T[days]
            cum_T10 = np.cumsum(T) / len(T)
            #cum_T10 = np.cumsum(T) / np.arange(1, len(T) + 1)

            cum_T10 = cum_T10[days]
            M_train.extend(M10)
            T_train.extend(cum_T10)
    if w0_test is None:
        print(f"{test_year}: skipped (no human cases)")
        continue

    # Convert to arrays for KDE
    M_train = np.asarray(M_train)
    T_train = np.asarray(T_train)
    # Fit a 2D KDE to the (M, T) pairs
    kde = gaussian_kde(np.vstack([M_train, T_train]))

    # Run compartmental model for test year using *average* temperature (excluding this year)
    avg_temp = average_temperature_excluding_year(test_year)
    n_days = len(avg_temp)
    C_ac_thisyear = Cac_interp(avg_temp)  # Get carrying capacity vector for each day
    Coematrix = build_coef(avg_temp)      # Get temperature-dependent model coefficients
                   
    x0  = np.array([3e5, 0, 0, 0, 1, 1, 100, 1, 1, 1e4, 1, 1, 1])  # Initial conditions for all compartments
    days = np.arange(n_days)
    ode_rhs = rhs_factory(C_ac_thisyear, Coematrix)
    sol = odeint(ode_rhs, x0, days)
    AM_model = ( sol[:, 6] + sol[:, 7] + sol[:, 8] )    # Modeled aquatic mosquito abundance
    
    cum_T_model = np.cumsum(avg_temp) / n_days  # Cumulative temperature, averaged for each day
    #cum_T_model = np.cumsum(avg_temp) / np.arange(1, len(avg_temp) + 1)
    



    # 2D grid for plotting the KDE PDF

    m_lo = min(0.1e8, AM_model.min()) * 0.2

 
    m_hi = max(M_train.max(), 1.65e8) * 1.2

    t_lo = min(T_train.min(), cum_T_model.min()) * 0.9
    t_hi = max(T_train.max(), cum_T_model.max()) * 1.1
    m_lin = np.linspace(m_lo, m_hi, 160)
    t_lin = np.linspace(t_lo, t_hi, 160)
    MG, TG = np.meshgrid(m_lin, t_lin)
    pdf = kde(np.vstack([MG.ravel(), TG.ravel()])).reshape(MG.shape)
    
    
    
    
    
    
    dx, dy = m_lin[1]-m_lin[0], t_lin[1]-t_lin[0]
    print("pdf min/max =", pdf.min(), pdf.max())
    print("approx integral over grid =", pdf.sum() * dx * dy)

    # Find the threshold for the 60% highest-probability region (alpha contour)
    dx, dy = m_lin[1]-m_lin[0], t_lin[1]-t_lin[0]
    flat   = pdf.ravel()
    order  = flat.argsort()[::-1]
    mass   = np.cumsum(flat[order]) * dx * dy
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    C_level = 0.5

    #max C_level = 0.999
    
    
    
    
    
    thr60 = flat[order[np.searchsorted(mass, C_level)]]

    # What were (M,T) on the first human case day in the test year?
    obs_day = w0_test * 7
    M_obs = AM_model[obs_day]
    T_obs = cum_T_model[obs_day]

    # ---- PLOT 1: Show KDE, alpha contour, and first case day ----
    plt.figure(figsize=(7,5))
    plt.contourf(MG, TG, pdf, levels=25, cmap='cividis')
    #plt.contour(MG, TG, pdf, levels=[thr60], colors='#E69F00', lw=2, ls='--')
    plt.contour(MG, TG, pdf, levels=[thr60], colors='red', lw=2, ls='-')
    plt.scatter(AM_model, cum_T_model, c='black', marker='.', s=8, edgecolor='black')#, label='$(M, \\bar T)$ at first case')
    plt.scatter([M_obs], [T_obs], c='white', marker='X', s=80, edgecolor='black', label='$(M, \\bar T)$ at first case')
    


    #cf = plt.contourf(MG, TG, pdf, levels=25, cmap='cividis')
    #cbar = plt.colorbar(cf, format="%.2e")  # scientific notation
    #cbar.set_label("KDE density")
        
    
    
    
    plt.colorbar(label='Density')
    plt.xlabel('Mosquito abundance (M)')
    plt.ylabel('Average total temperature (°C·day)')
    plt.title(f'Orange County {test_year}-confidence level={C_level}')
    plt.legend(loc='upper right', fontsize=12)
    plt.tight_layout()
    plt.xlim(0.4e8,1.9e8)
    plt.ylim(6,14)
    plt.savefig(f'interval_estimation_retro_PDF_{C_level}_{test_year}_or.jpg', format='jpg', bbox_inches='tight')

    plt.show()

    # ---- PLOT 2: Timeline of days flagged as high risk by PDF ----
    # For every day, check if model trajectory is inside the "high risk" contour
    kde_vals = kde(np.vstack([AM_model, cum_T_model]))
    inside = kde_vals >= thr60
    red_days = np.where(inside)[0]
    red_days = red_days[ red_days >= min_index ] 
    #print(red_days)
    #select_block_or_all(red_days, obs_day)
    #print("second tr",select_block_or_all(red_days, obs_day)[1])
    if select_block_or_all(red_days, obs_day)[1]==1:
        red_days=select_block_or_all(red_days, obs_day)[0]
        #print(red_days)
        inside[:]=False
        inside[red_days]=True
        #print("R",red_days)
    
        #N_blocks=select_block_or_all(red_days, obs_day)[1]
    
        if red_days.size:
    
            if (len(red_days)) > 1:
                #print(obs_day)
                #print(f"{test_year} → red days: {red_days[0]}, {red_days[-1]}")
                if obs_day >= red_days[0] and obs_day <= red_days[-1]:

                    #print("length of red",len(inside[red_days]))
                    information.append([test_year,len(inside[red_days]),np.round(1/len(inside[red_days]),decimals=3)])
    
                    u += 1
                    plt.figure(figsize=(7,5))
                    for d, flag in enumerate(inside):
                        col = '#E69F00' if flag else '#56B4E9'  # Orange if "in high risk", sky blue otherwise
                        plt.axvline(d, 0, 1, color=col, lw=2)
                    plt.axvline(obs_day, 0, 1, color='black', lw=2, label='First human-case day')
                    plt.xlim(0, len(inside)-1)
                    plt.yticks([])
                    plt.xlabel('Day of year')
                    plt.title(f'Risk timeline – {test_year}-confidence level={C_level}-Orange County')
                    plt.legend(loc='upper left', fontsize=12)
                    plt.tight_layout()
                 
                    plt.savefig(f'interval_estimation_retro_{C_level}_{test_year}_or.jpg', format='jpg', bbox_inches='tight')

                    plt.show()
    else:
        information.append([test_year,len(inside[red_days]),0])

                        
                    
        plt.figure(figsize=(7,5))
        for d, flag in enumerate(inside):
            col = '#E69F00' if flag else '#56B4E9'  # Orange if "in high risk", sky blue otherwise
            plt.axvline(d, 0, 1, color=col, lw=2)
        plt.axvline(obs_day, 0, 1, color='black', lw=2, label='First human-case day')
        plt.xlim(0, len(inside)-1)
        plt.yticks([])
        plt.xlabel('Day of year')
        plt.title(f'Risk timeline – {test_year}-confidence level={C_level}-Orange County')
        plt.legend(loc='upper left', fontsize=12)
        plt.tight_layout()
    
        plt.savefig(f'interval_estimation_retro_{C_level}_{test_year}_or.jpg', format='jpg', bbox_inches='tight')
        plt.show()

        

# ---- PLOT: Lead time summary ----

print("confidence level $alpha$=",C_level)
print("ability-$beta$",u/len(YEARS))
print("year,interval length, accuracy")
#print(np.array(information,dtype=float))

arr = np.array(information, dtype=float)
print(np.array2string(
    arr,
    formatter={'float_kind': lambda x: f'{x:.3f}'},  # set your decimals
    max_line_width=200
))



def plot_inverse_length(info):
    """
    Parameters
    ----------
    info : list[list]  –  rows of [year, length, error (=1/length)]
                         order can be arbitrary
    """
    if not info:
        raise ValueError("`info` is empty")

    # --- sort by year --------------------------------------------------
    info_sorted = sorted(info, key=lambda row: row[0])

    years  = [int(row[0]) for row in info_sorted]   # ensure plain ints
    errors = [np.round(row[2],decimals=2)          for row in info_sorted]

    # --- plot ----------------------------------------------------------
    plt.figure(figsize=(8, 4))
    plt.plot(years, errors, 'o-', markerfacecolor='white',
             markeredgecolor='black', lw=1.5)
    plt.xticks(years, rotation=45)                  # integer ticks
    plt.xlabel('Year')
    plt.ylabel('Error')
    plt.title('Inverse-length error by season')
    plt.tight_layout()
    plt.show()

# ----------------------------------------------------------------------
# call it
plot_inverse_length(information)

from ast import literal_eval
from typing import Union, Iterable, List, Dict

def compute_avgs(information: Union[str, Iterable[Iterable[float]]]) -> Dict[str, float]:
    """
    information: either a string like '[[a,b,c], ...]' or an iterable of [a,b,c].
    Returns:
        {
          "avg_L_alpha": average of all b's,
          "avg_a_alpha": average of all c's
        }
    """
    # Parse if given as a string
    data = literal_eval(information) if isinstance(information, str) else list(information)
    if not data:
        raise ValueError("information is empty")

    b_vals = [row[1] for row in data]
    c_vals = [row[2] for row in data]

    avg_L_alpha = sum(b_vals) / len(b_vals)
    avg_a_alpha = sum(c_vals) / len(c_vals)

    return {"avg_L_alpha": avg_L_alpha, "avg_a_alpha": avg_a_alpha}


# --- Example with your data ---

results = compute_avgs(information)
# print(results)  # {'avg_L_alpha': 26.736842105263158, 'avg_a_alpha': 0.027038076876493827}



