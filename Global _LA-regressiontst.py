# -*- coding: utf-8 -*-
"""
Created on Mon Jun 16 09:30:30 2025

@author: shosseini
"""

# -*- coding: utf-8 -*-
"""
Created on Fri Jun 13 16:10:09 2025

@author: shosseini
"""


import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from scipy.stats import gaussian_kde
from scipy.ndimage import gaussian_filter1d

# Y = [2005,
#     2006, 2007, 2008, 2009,2010,
#       2011, 2012, 2013,
#     2015, 2016, 2017, 2018,
#     2019, 2020, 2021, 2022,
#     2023,2024
# ]
Y = [y for y in range(1990, 2025) ]#if y not in (2010, 2011)]

point_estimation=[]
Global=[]
for exy in Y:

    print("-----------------------------------------excluded year",exy)

    YEARS = [
        2006, 2007, 2008, 2009,2010,
        2011, 2012, 2013,2014,
        2015, 2016, 2017, 2018,
        2019, 2020, 2021,2022,
        2023,2024
    ]
    
    # Calculate total number of cases across all years
    total_cases = 0
    for year in YEARS:
        try:
            

            weekly_cases = np.loadtxt(f'LA{year}.txt', dtype=int)
            
            
            if weekly_cases.size < 52:
                weekly_cases = np.pad(weekly_cases, (0, 52 - weekly_cases.size), constant_values=0)
            total_cases += weekly_cases.sum()
        except Exception as e:
            print(f"Skipping year {year} for total cases count due to error: {e}")
    
    #print(f"Total number of human cases across all years: {total_cases}")
    
    # Load climate data once
    climate_df = pd.read_csv('LA1.csv')

    
    
    def first_nonzero_index(year):
        

        lst=list(np.loadtxt(f'LA{year}.txt', dtype=int))
        
        
        
        for i in range(len(lst)):
            #print(lst[i])
            if lst[i] > 0:
                a=7*i
                break
        return a
        
    
    
 
    
    

    
    
    
    
    

    
    
    
    def temperature_series(year):
        temps = climate_df.loc[climate_df['Year'] == year, 'T'].to_numpy()
        if temps.size == 0:
            raise ValueError(f"No temperature data for Year={year}")
        return temps
    
    def build_weighted_MT_excluding_year(years):
        #MSE_list=[]
        """
        For each week with n cases, collect the 7 days from 7 to 0 days prior to the start of that week,
        and add them to the data n times. Also returns the maximum repetition count of any (M, T) pair.
        """
        M_vals = []
        T_vals = []
        for yr in years:
            

            weekly_cases = np.loadtxt(f'LA{yr}.txt', dtype=int)
            
            
            if weekly_cases.size < 52:
                weekly_cases = np.pad(weekly_cases, (0, 52 - weekly_cases.size), constant_values=0)
                
            M_daily = np.loadtxt(f'MLA_{yr}.txt')

            
            
            T_daily = temperature_series(yr)
            length = min(len(M_daily), len(T_daily))
            M_daily = M_daily[:length]
            T_daily = np.cumsum(T_daily[:length]) / length
            for w, cases in enumerate(weekly_cases):
                if cases == 0:
                    continue
                # 7 days: from 7 days before to the day week start
                start = max(0, w * 7 - 14)
                end = max(0, w * 7-7)  # exclusive, so stops at w*7
                M_window = M_daily[start:end]
                T_window = T_daily[start:end]
                for _ in range(cases):
                    M_vals.extend(M_window)
                    T_vals.extend(T_window)

    
        M_arr = np.array(M_vals)
        T_arr = np.array(T_vals)
    
        # Calculate max repetition of any (M, T) pair
        if len(M_arr) > 0:
            # Use np.round to deal with float precision if needed; here using 2 decimals
            pairs = np.column_stack((np.round(M_arr, 2), np.round(T_arr, 3)))
            unique_pairs, counts = np.unique(pairs, axis=0, return_counts=True)
            max_repetition = counts.max()
            #print(f"Maximum repetition of a single (M, T) pair: {max_repetition}")
        else:
            max_repetition = 0
            print("No data collected.")
    
        return M_arr, T_arr, max_repetition
    
    
    
    def plot_excluded_year_kde_curve(years):
        M_train, T_train = build_weighted_MT_excluding_year(years)[0:2]
        if len(M_train) == 0:
            print(f"No training data available excluding year {excluded_year}.")
            return
    
        # KDE PDF (do not scale)
        data_train = np.vstack([M_train, T_train])
        kde = gaussian_kde(data_train)#, bw_method='scott')
    
        x_min, x_max = 0, M_train.max() + 0.5e8
        y_min, y_max = 3, 20
        xbins, ybins = 100,100
        x_grid = np.linspace(x_min, x_max, xbins)
        y_grid = np.linspace(y_min, y_max, ybins)
        X, Y = np.meshgrid(x_grid, y_grid)
    
        # This is the PDF, which integrates to 1 (approximately, numerically)
        kde_pdf = kde(np.vstack([X.ravel(), Y.ravel()])).reshape(X.shape)
        correctionvalue=build_weighted_MT_excluding_year(years)[2]/np.max(kde_pdf)
        kde_pdf = correctionvalue*kde(np.vstack([X.ravel(), Y.ravel()])).reshape(X.shape)
        
        # Print peak value of the PDF
        print(f"Peak value of the PDF (maximum density): {np.max(kde_pdf):.4g}")
    
        # Print numerical integral (should be close to 1)
        print("Numerical integral (should be close to 1):", 
              np.sum(kde_pdf) * (x_grid[1]-x_grid[0]) * (y_grid[1]-y_grid[0]))
    
        cmap = plt.get_cmap('cividis')
        plt.figure(figsize=(14, 8))
        cf = plt.contourf(X, Y, kde_pdf, levels=25, cmap=cmap)
        plt.colorbar(cf, label='Probability Density (PDF)')
        plt.xlabel('Mosquito abundance (M)')
        plt.ylabel('Mean cumulative temperature (°C·day)')
        plt.title(f'Year {excluded_year}: KDE Probability Density Function (PDF)')
        plt.tight_layout()
        plt.show()
    
    
    
    
    
    
    
    
        
        
        
    excluded_year = exy
    # plot_excluded_year_kde_curve(YEARS)
    
    def plot_excluded_year_kde_curve_and_daily_risk(years,excluded_year,C_predict_file):
        MSE_list=[]
        Gl=0
        """
        Long-term cross-validation plot for a leave-one-year-out experiment.
        ▸ Top: scaled KDE surface + (M,T) trajectory of the excluded season
        ▸ Bottom: daily-risk timeline coloured with the SAME scale
        """
    
        import numpy as np
        import matplotlib.pyplot as plt
        import matplotlib.colors as mcolors
        from scipy.stats import gaussian_kde
        from scipy.integrate import odeint
        from scipy.interpolate import interp1d
    
    
        
    

    
    
    
        # ──────────────────────────────────────────────────────────────
        # Build temperature-dependent coefficient matrix
        # (unchanged from your original)
        # ──────────────────────────────────────────────────────────────
        def build_coef(T_series):
            def fwrap(fun):
                return np.vectorize(lambda d: fun(T_series[int(d)]))
            gamma = fwrap(lambda T: 0 if (T <= 4.3 or T >= 39.9)
                                    else T*4.12e-5*(T-4.3)*np.sqrt(39.9-T))
            pLST  = fwrap(lambda T: 0 if (T <= 5.9 or T >= 43.1)
                                    else -2.12e-3*(T-5.9)*(T-43.1))
            muM   = fwrap(lambda T: 1/1e-16 if T > 41.3
                                    else (1/45 if T < 14
                                          else 1/(0.99*(-1.69*T+69.6))))
            BRT   = fwrap(lambda T: 0 if (T <= 2.3 or T >= 32)
                                    else T*1.67e-4*(T-2.3)*np.sqrt(32-T))
            alpha = fwrap(lambda T: 0 if (T <= 5.3 or T >= 38.9)
                                    else -5.98e-1*(T-5.3)*(T-38.9)*BRT(int(T)))
            bET   = fwrap(lambda T: 0 if (T <= 3.2 or T >= 42.6)
                                    else -2.11e-3*(T-3.2)*(T-42.6))
            PDRT  = fwrap(lambda T: 0 if (T <= 11.2 or T >= 44.7)
                                    else T*6.57e-5*(T-11.2)*np.sqrt(44.7-T))
            VCT   = fwrap(lambda T: 0 if (T <= 11.3 or T >= 41.9)
                                    else -2.94e-3*(T-11.3)*(T-41.9))
            HR     = lambda d: 1.0
            muegg  = lambda d: 0.0
            ADR    = lambda d: pLST(d)*gamma(d)
            muAQ   = lambda d: gamma(d)*(1-pLST(d))
            funcs = (gamma, pLST, muM, BRT, alpha, bET,
                      PDRT, VCT, HR, muegg, ADR, muAQ)
    
            M = np.zeros((12, 365))
            for day in range(len(T_series)):
                M[:, day] = [f(day) if not callable(f) else f(day) for f in funcs]
            return M
    
        # ──────────────────────────────────────────────────────────────
        # ODE RHS factory (unchanged)
        # ──────────────────────────────────────────────────────────────
        def rhs_factory(C_ac_vec, Coe, control=0, controlday=0):
            def BMPPT(x, t):
                day = int(t) % 365
                C_ac = C_ac_vec[int(t)]
                muB = 0.15/365; muhum = 0.01; incH = 1/6; infH = 0.14
                muH = 1/(77*365); incB = 1/3; recC = 1/6; muwnvB = 0.9
                SH, EH, IH, RH, EM, AM, MS, ME, MI, BS, BE, BI, BR = x
                TH = SH+EH+IH+RH; D = BS+BE+BI+BR; eps = 1e-8
                ph = (Coe[3, day] * Coe[7, day]) / (D + TH + eps)
                # derivatives …
                degg = Coe[4,day]*(MS+ME+MI) - Coe[5,day]*EM \
                        - Coe[8,day]*EM - Coe[9,day]*EM
                dacu = Coe[5,day]*EM*max(0, 1-AM/(C_ac+eps)) \
                        - Coe[0,day]*Coe[1,day]*AM \
                        - Coe[10,day]*AM - Coe[11,day]*AM
                return np.array([
                    muH*TH - ph*MI*SH      - muH*SH,                # dSH
                    ph*MI*SH - incH*EH     - muH*EH,                # dEH
                    incH*EH - infH*IH - muhum*IH - muH*IH,          # dIH
                    infH*IH - muH*RH,                               # dRH
                    degg, dacu,
                    Coe[10,day]*AM - ph*BI*MS - Coe[2,day]*MS,      # dMS
                    ph*BI*MS - Coe[6,day]*ME - Coe[2,day]*ME,       # dME
                    Coe[6,day]*ME - Coe[2,day]*MI,                  # dMI
                    muB*D - ph*Coe[7,day]*MI*BS,                    # dBS
                    ph*Coe[7,day]*MI*BS - incB*BE - muB*BE,         # dBE
                    incB*BE - muwnvB*BI - muB*BI - recC*BI,         # dBI
                    recC*BI - muB*BR                                # dBR
                ])
            return BMPPT
    
        # ──────────────────────────────────────────────────────────────
        # 1.  Build KDE (training data excluding test year)
        # ──────────────────────────────────────────────────────────────
        minindex=1#first_nonzero_index(excluded_year)
        #print("min",minindex)
        
        M_train, T_train, max_rep = build_weighted_MT_excluding_year(years)
        if len(M_train) == 0:
            print(f"No training data available excluding year {excluded_year}.")
            return
    
        kde = gaussian_kde(np.vstack([M_train, T_train]), bw_method='scott')
    
        x_min, x_max = 0, M_train.max() + 0.5e8
        y_min, y_max = 3, 20
        xbins = ybins = 100
        x_grid = np.linspace(x_min, x_max, xbins)
        y_grid = np.linspace(y_min, y_max, ybins)
        X, Y = np.meshgrid(x_grid, y_grid)
    
        # Raw PDF
        pdf_vals = kde(np.vstack([X.ravel(), Y.ravel()]))
        # Scale so that max(pdf_vals) == max repetition count
        correction_value = max_rep / np.max(pdf_vals)
        counts_kde = (correction_value * pdf_vals).reshape(X.shape)  # “expected cases”
    
        norm = mcolors.Normalize(vmin=0, vmax=counts_kde.max())
        cmap = plt.get_cmap('cividis')
    
        # ──────────────────────────────────────────────────────────────
        # 2.  Prepare ODE simulation for the excluded year
        # ──────────────────────────────────────────────────────────────
        avg_temp = temperature_series(excluded_year)[0:365]
        Coe      = build_coef(avg_temp)
    
        # Carrying-capacity function
        C_pred   = np.loadtxt(C_predict_file, delimiter=',', skiprows=1)
        T_vals   = C_pred[:, 0];    C_vals = C_pred[:, 1]
        C_interp = interp1d(T_vals, C_vals, kind='linear',
                            bounds_error=False, fill_value='extrapolate')
        C_ac_daily = C_interp(avg_temp[:364])
    
        rhs = rhs_factory(C_ac_daily, Coe)
        x0  = np.array([3e5, 0, 0, 0, 1, 1, 100, 1, 1, 1e4, 1, 1, 1])
        sol = odeint(rhs, x0, np.arange(364))
        M_y =list( sol[:, 6] + sol[:, 7] + sol[:, 8] )         # total adult mosquitoes)
        cum_T_y = np.cumsum(avg_temp[:364]) / len(avg_temp[:364])
    
        # ──────────────────────────────────────────────────────────────
        # 3.  Figure with two stacked panels
        # ──────────────────────────────────────────────────────────────
        fig, (ax1, ax2) = plt.subplots(2, 1, figsize=(14, 8),
                                        gridspec_kw={'height_ratios':[1, 1]})
    
        # ── top panel ────────────────────────────────────────────────
        cf = ax1.contourf(X, Y, counts_kde, levels=25, cmap=cmap, norm=norm)
        fig.colorbar(cf, ax=ax1).set_label('Expected case count', fontsize=14)
    
        ax1.plot(M_y, cum_T_y, c='cyan', lw=3,
                  label=f'{excluded_year}  (M, T)')
        ax1.set_xlabel('Mosquito abundance (M)', fontsize=14)
        ax1.set_ylabel('T', fontsize=14)
        ax1.set_title(f'Retrospective forecast for {excluded_year}', fontsize=15)
        ax1.legend()
        ax1.set_ylim(8.5, y_max)
    
        # ── bottom panel : daily-risk stripes ────────────────────────
        daily_intensity = []
        for m_val, t_val in zip(M_y, cum_T_y):
            if (m_val < x_min) or (m_val > x_max) or (t_val < y_min) or (t_val > y_max) or M_y.index(m_val)<minindex:
                daily_intensity.append(0)                     # outside domain
            else:
                
                dens = kde(np.array([[m_val], [t_val]]))[0]
                daily_intensity.append(dens * correction_value)
        daily_intensity = np.array(daily_intensity)
    
    
        for d, inten in enumerate(daily_intensity):
            if inten>=1:
                Gl=Gl+1
            #col = '#56B4E9' if inten < 0 else cmap(norm(inten))
            #print(d)
            col=cmap(norm(inten))
            ax2.axvline(d, 0, 1, color=col, lw=4)
            
            
            
        Global.append(Gl)

    
        ax2.set_xlim(0, len(M_y))
        ax2.set_ylim(0, 1)
        ax2.set_yticks([])
        ax2.set_xlabel('Day of year', fontsize=14)
    
        sm = plt.cm.ScalarMappable(cmap=cmap, norm=norm)
        sm.set_array([])
        fig.colorbar(sm, ax=ax2, orientation='vertical', pad=0.02)\
            .set_label('Expected case count', fontsize=14)
    
        plt.tight_layout()
        # plt.savefig(f'Retrospective_Severity_map_{excluded_year}_orange-1year.jpg', format='jpg', bbox_inches='tight')
        # plt.savefig(f'Retrospective_Severity_map_{excluded_year}_orange-1year.pdf', format='pdf', bbox_inches='tight')
        plt.show()
        #print(MSE_list)
        import numpy as np
    

    
    
    # Assume all the imports, climate_df loading, and helper functions are defined above
    
    excluded_year = excluded_year # The year you want to exclude and predict
    
    # Call the plotting function
    
    
    
    
    

    MSE_list=plot_excluded_year_kde_curve_and_daily_risk(YEARS, excluded_year, C_predict_file='C_predict_60_LA.txt')
    
    
    
    
    
    
plt.plot(Global)
plt.show()    
    





from scipy import stats
import pymannkendall as mk  # Ensure this is installed: pip install pymannkendall

def fit_and_plot_linear_trend(ts_values, time_index=None):
    """
    Fits a linear regression to the given time series and plots the result.
    Also performs additional statistical tests to assess the presence of a trend.
    
    Parameters:
    - ts_values: list or array of time series values
    - time_index: optional list or array of time points (e.g., years). If None, uses 0, 1, 2, ...
    
    Returns:
    - slope, intercept, r_squared, p_value
    """
    y = np.array(ts_values)
    
    if time_index is None:
        x = np.arange(len(y))
    else:
        x = np.array(time_index)

    # Linear regression
    slope, intercept, r_value, p_value, std_err = stats.linregress(x, y)
    regression_line = intercept + slope * x

    # Additional test: Spearman's rank correlation
    spearman_r, spearman_p = stats.spearmanr(x, y)

    # Additional test: Mann-Kendall trend test
    mk_result = mk.original_test(y)

    # Print regression results
    print(f"Linear Regression:")
    print(f"  Slope: {slope:.4f}")
    print(f"  Intercept: {intercept:.4f}")
    print(f"  R-squared: {r_value**2:.4f}")
    print(f"  P-value: {p_value:.4e}")
    print()

    print(f"Spearman's Rank Correlation:")
    print(f"  Spearman's rho: {spearman_r:.4f}")
    print(f"  P-value: {spearman_p:.4e}")
    print()

    print(f"Mann-Kendall Trend Test:")
    print(f"  Trend: {mk_result.trend}")
    print(f"  P-value: {mk_result.p:.4e}")
    print(f"  Slope: {mk_result.slope:.4f}")
    print()


    
    
    plt.figure(figsize=(8, 4))
    plt.plot(x, y, label="${N_i}$ time series", marker='o')
    plt.plot(x, regression_line, label="Linear Trend", linestyle='--')
    plt.title("Time series of the number of risky days with linear trend for Los Angeles County, California")
    plt.xlabel("Time" if time_index is None else "Year")
    plt.ylabel("Estimated number of risky days")
    plt.xticks(ticks=x, labels=x, rotation=45)  # <-- This is the only added line
    plt.legend()
    plt.grid(False)
    plt.tight_layout()
    plt.savefig('trend_Los Angeles.pdf', format="pdf", bbox_inches="tight")
    plt.savefig('trend_Los Angeles.jpg', format="jpg", bbox_inches="tight")
    plt.show()

    return slope, intercept, r_value**2, p_value











fit_and_plot_linear_trend(Global, Y)   
