# -*- coding: utf-8 -*-
"""
Created on Fri Jun 13 16:10:09 2025

@author: shosseini
"""


import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import matplotlib.ticker as mtick
import math
import os

Y = [
    2006, 2007
    ,2008, 2009,2010,
      2011, 2012, 2013,
    2015, 2016, 2017, 2018,
    2019, 2020, 2021, 2022,
    
    2023
    ,
    2024
]

point_estimation=[]
em=0
tm=0
T_log_scores=[]
comparison_MSE=[]
Comaoprison_realMSE=[]
for exy in Y:
    print("-----------------------------------------excluded year",exy)

    YEARS = [
        2006, 2007, 2008, 2009,2010,
         2011, 2012, 2013,
        2015, 2016, 2017, 2018,
        2019, 2020, 2021, 2022,
        2023
        ,2024
    ]
    
    # Calculate total number of cases across all years
    total_cases = 0
    for year in YEARS:
        try:
            weekly_cases = np.loadtxt(f'orangecases{year}.txt', dtype=int)
            if weekly_cases.size < 52:
                weekly_cases = np.pad(weekly_cases, (0, 52 - weekly_cases.size), constant_values=0)
            total_cases += weekly_cases.sum()
        except Exception as e:
            print(f"Skipping year {year} for total cases count due to error: {e}")
    
    #print(f"Total number of human cases across all years: {total_cases}")
    
    # Load climate data once
    climate_df = pd.read_csv('orange.csv')
    
    
    def first_nonzero_index(year):
        lst=list(np.loadtxt(f'orangecases{year}.txt', dtype=int))
        for i in range(len(lst)):
            #print(lst[i])
            if lst[i] > 0:
                a=7*i
                break
        return a
        
    
    
    from scipy.stats import gaussian_kde
    
    

    
    
    
    
    
    def mean_r_excluding(excluded_year, folder_path="."):
        """
        Compute the mean of all weekly values in OR{year}.txt files,
        excluding 2010, 2022, and the specified excluded_year.
        
        Parameters:
            excluded_year (int): Year to exclude from averaging.
            folder_path (str): Directory containing OR{year}.txt files.
        
        Returns:
            float: Mean of all values across all eligible years.
        """
        valid_years = [y for y in range(2006, 2025) if y not in {2014, excluded_year}]
        all_values = []
    
        for year in valid_years:
            file_path = os.path.join(folder_path, f"orangecases{year}.txt")
            if os.path.exists(file_path):
                try:
                    data = np.loadtxt(file_path, dtype=int)
                    all_values.extend(data)
                except Exception as e:
                    print(f"Warning: could not read OR{year}.txt: {e}")
            else:
                print(f"Warning: file OR{year}.txt not found")
    
        if not all_values:
            raise ValueError("No valid data found across available years.")
    
        return np.mean(all_values)
    
    
    
    # Example usage (disabled actual execution for safety):
    # mean_array = mean_rs_files_excluding(2015, folder_path="path/to/files")
    # mean_array
    
    
    
    def temperature_series(year):
        temps = climate_df.loc[climate_df['Year'] == year, 'T'].to_numpy()
        if temps.size == 0:
            raise ValueError(f"No temperature data for Year={year}")
        return temps
    
    def build_weighted_MT_excluding_year(years, excluded_year):
        #MSE_list=[]
        """
        For each week with n cases, collect the 7 days from 13 to 7 days prior to the start of that week,
        and add them to the data n times. Also returns the maximum repetition count of any (M, T) pair.
        """
        M_vals = []
        T_vals = []
        for yr in years:
            if yr == excluded_year:
                continue
            try:
                weekly_cases = np.loadtxt(f'orangecases{yr}.txt', dtype=int)
                if weekly_cases.size < 52:
                    weekly_cases = np.pad(weekly_cases, (0, 52 - weekly_cases.size), constant_values=0)
                M_daily = np.loadtxt(f'MOR_{yr}.txt')
                T_daily = temperature_series(yr)
                length = min(len(M_daily), len(T_daily))
                M_daily = M_daily[:length]
                T_daily = np.cumsum(T_daily[:length]) / length
                for w, cases in enumerate(weekly_cases):
                    if cases == 0:
                        continue
                    # 7 days: from 7 days before to the day week start
                    start = max(0, w * 7 - 14)
                    end = max(0, w * 7 -7)  # exclusive, so stops at w*7
                    M_window = M_daily[start:end]
                    T_window = T_daily[start:end]
                    for _ in range(cases):
                        M_vals.extend(M_window)
                        T_vals.extend(T_window)
            except Exception as e:
                print(f"Skipping year {yr} due to error: {e}")
    
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
    
    
    
    def plot_excluded_year_kde_curve(years, excluded_year):
        M_train, T_train = build_weighted_MT_excluding_year(years, excluded_year)[0:2]
        if len(M_train) == 0:
            print(f"No training data available excluding year {excluded_year}.")
            return
    
        # KDE PDF (do not scale)
        data_train = np.vstack([M_train, T_train])
        kde = gaussian_kde(data_train, bw_method='scott')
    
        x_min, x_max = 0, M_train.max() + 0.5e8
        y_min, y_max = 3, 20
        xbins, ybins = 100,100
        x_grid = np.linspace(x_min, x_max, xbins)
        y_grid = np.linspace(y_min, y_max, ybins)
        X, Y = np.meshgrid(x_grid, y_grid)
    
        # This is the PDF, which integrates to 1 (approximately, numerically)
        kde_pdf = kde(np.vstack([X.ravel(), Y.ravel()])).reshape(X.shape)
        correctionvalue=build_weighted_MT_excluding_year(years, excluded_year)[2]/np.max(kde_pdf)
        kde_pdf = correctionvalue*kde(np.vstack([X.ravel(), Y.ravel()])).reshape(X.shape)
        
        # Print peak value of the PDF
        #print(f"Peak value of the PDF (maximum density): {np.max(kde_pdf):.4g}")
    
        # Print numerical integral (should be close to 1)
        # print("Numerical integral (should be close to 1):", 
        #       np.sum(kde_pdf) * (x_grid[1]-x_grid[0]) * (y_grid[1]-y_grid[0]))
    
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
    plot_excluded_year_kde_curve(YEARS, excluded_year)
    
    def plot_excluded_year_kde_curve_and_daily_risk(years,excluded_year,C_predict_file='C_predict_75.txt'):
        MSE_list=[]
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
    
    
        
    
        # def average_temperature(years_list,excluded):
        #     """
        #     Return the daily mean-temperature series obtained from the five seasons
        #     immediately preceding *excluded*. Requires climate_df to contain columns 'Year' and 'T'.
        #     """
        #     prev_years = [excluded - i for i in range(1, 6)]
        #     temps = []
        #     for y in prev_years:
        #         t = climate_df.loc[climate_df['Year'] == y, 'T'].to_numpy()
        #         if t.size:
        #             temps.append(t[:364])
        #         else:
        #             print(f"Warning: no temperature data for {y}")
        #     if len(temps) == 0:
        #         raise RuntimeError("No temperature data available for the preceding years")
        #     return np.mean(temps, axis=0)
        
        
        
        
        def average_temperature(years_list,excluded):
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
            head = full_year[:220]

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
            tail = clim_mean[220:]

            # ----- 3. concatenate head + tail ----------------------------------------
            return np.concatenate([head, tail])
    
    
    
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
        minindex=first_nonzero_index(excluded_year)
        #print("min",minindex)
        
        M_train, T_train, max_rep = build_weighted_MT_excluding_year(years, excluded_year)
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
        avg_temp = average_temperature(years, excluded_year)
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
        ax1.set_title(f'Cross-validation – {excluded_year} (long-term)', fontsize=15)
        ax1.legend()
        ax1.set_ylim(8.5, y_max)
    
        # ── bottom panel : daily-risk stripes ────────────────────────
        daily_intensity = []
        for m_val, t_val in zip(M_y, cum_T_y):
            if (m_val < x_min) or (m_val > x_max) or (t_val < y_min) or (t_val > y_max) or M_y.index(m_val)<minindex+5:
                daily_intensity.append(-1)                     # outside domain
            else:
                
                dens = kde(np.array([[m_val], [t_val]]))[0]
                daily_intensity.append(dens * correction_value)
        daily_intensity = np.array(daily_intensity)
    
    
        for d, inten in enumerate(daily_intensity):
            col = '#56B4E9' if inten < 0 else cmap(norm(inten))
            #print(d)
            ax2.axvline(d, 0, 1, color=col, lw=4)
    
        # Overlay weekly case counts for the excluded year
        try:
            weekly_cases = np.loadtxt(f'orangecases{excluded_year}.txt', dtype=int)
            if weekly_cases.size < 52:
                weekly_cases = np.pad(weekly_cases,
                                      (0, 52-weekly_cases.size), constant_values=0)
            for wk, n in enumerate(weekly_cases):
                if n>0 and wk*7>minindex:
                    
                    
                    
                    #print(wk*7,n,daily_intensity[wk*7])
                    
                    
                    
                    MSE_list.append((wk*7,n,daily_intensity[wk*7]))
                    point_estimation.append([n,daily_intensity[wk*7]])
                    
                    
                    
                if n > 0:
                    day_pos = wk*7
                    if day_pos>minindex:
                        ax2.axvline(day_pos, 0, 1, color='white', lw=2)
                        ax2.text(day_pos+4, 0.5, str(n), color='white',fontsize=14, ha='center', va='center', rotation=90,fontweight='bold')
                    #print(wk*7,weekly_cases[wk*7],daily_intensity[wk*7])
        except Exception as e:
            print(f'Weekly case file missing for {excluded_year}: {e}')
    
        ax2.set_xlim(0, len(M_y))
        ax2.set_ylim(0, 1)
        ax2.set_yticks([])
        ax2.set_xlabel('Day of year', fontsize=14)
    
        sm = plt.cm.ScalarMappable(cmap=cmap, norm=norm)
        sm.set_array([])
        fig.colorbar(sm, ax=ax2, orientation='vertical', pad=0.02)\
           .set_label('Expected case count', fontsize=14)
    
        plt.tight_layout()
        plt.savefig(f'Severity_map_{excluded_year}_orange-1year.jpg', format='jpg', bbox_inches='tight')
        plt.savefig(f'Severity_map_{excluded_year}_orange-1year.pdf', format='pdf', bbox_inches='tight')
        plt.show()
        #print(MSE_list)
        import numpy as np
    
        def rmse_day_obs_pred(triples):
            """
            triples : list of (day, observation, prediction) tuples
        
            Calculates RMSE using only the observation (index 1) and
            prediction (index 2) values, ignoring the day stamp
            (index 0).
            """
            if not triples:
                raise ValueError("Input list is empty")
        
            # extract observations and predictions
            obs  = np.array([t[1] for t in triples], dtype=float)
            pred = np.array([t[2] for t in triples], dtype=float)
         
            prednul=np.array([mean_r_excluding(excluded_year, folder_path=".") for t in triples], dtype=float)
    
        
            # root-mean-square error
            return float(np.sqrt(np.mean((obs - pred) ** 2))/np.mean(pred)),float(np.sqrt(np.mean((obs - prednul) ** 2))/np.mean(prednul)),float(np.sqrt(np.mean((obs - pred) ** 2)))
    
        
        
        # ---- example usage -------------------------------------------------
        RMSE = rmse_day_obs_pred(MSE_list)[0]
        print(f"RMSE = {RMSE:.4f}")
        RMSENUL=rmse_day_obs_pred(MSE_list)[1]
        print("RMSE_null",RMSENUL)
        comparison_MSE.append([RMSE,RMSENUL])
        Comaoprison_realMSE.append(rmse_day_obs_pred(MSE_list)[2])
        
        return(MSE_list,comparison_MSE, Comaoprison_realMSE)
    
    
    
    
    # Assume all the imports, climate_df loading, and helper functions are defined above
    
    # excluded_year = excluded_year # The year you want to exclude and predict
    
    # Call the plotting function
    MSE_list=plot_excluded_year_kde_curve_and_daily_risk(YEARS, excluded_year, C_predict_file='C_predict.txt')[0]
    
    
    
    
    
    
    
    # --- bin definitions -------------------------------------------------
    bins = [
        (0, 0), (1, 5), (6, 10), (11, 15), (16, 20),
        (21, 25), (26, 30), (31, 35), (36, 40), (41, 45),
        (46, 50), (51, 100), (101, 150), (151, 200), (201, 250)
    ]
    
    def bin_index_for_value(val):
        for idx, (lo, hi) in enumerate(bins):
            if lo <= val <= hi:
                return idx
        return None
    
    def poisson_pmf(k, lam):
        k_arr = np.asarray(k, dtype=float)
        lgamma_vals = np.array([math.lgamma(int(ki) + 1) for ki in k_arr])
        return np.exp(k_arr * np.log(lam) - lam - lgamma_vals)
    
    def bin_probability(lam, idx):
        lo, hi = bins[idx]
        ks = np.arange(lo, hi + 1)
        return poisson_pmf(ks, lam).sum()
    
    log_scores = []
    
    for week_num, (day, obs, lam) in enumerate(MSE_list, 1):
        # lam=np.round(lam)
        idx = bin_index_for_value(obs)
        if idx is None:
            score = -10.0
        else:
            P = bin_probability(lam, idx)
            score = np.log(P) if P > 0 else -10.0
            score = max(score, -10.0)
        log_scores.append(score)
    
        # ----------------- plotting --------------------------------------
        k_vals  = np.arange(0, 20)
        pmf_vals = poisson_pmf(k_vals, lam)
    
        fig, ax = plt.subplots(figsize=(14, 8))
        ax.plot(k_vals, pmf_vals, '--', linewidth=4)
    
        # --- NEW: shade the bin containing the observation ---------------
        if idx is not None:
            lo, hi = bins[idx]
            mask = (k_vals >= lo) & (k_vals <= hi)
            ax.fill_between(k_vals[mask], 0, pmf_vals[mask],
                            color='green', alpha=0.35, label=f'-log(Observed‐bin area)={score:.2f}', )
    
        ax.set_xlabel("Number of cases", fontsize=16)
        ax.set_ylabel("P(k | λ)", fontsize=16)
        ax.set_title(
            f"Week {day/7} ({excluded_year}): Predicted Poisson Distribution (λ = {lam:.2f}) "
            "for Orange County, California",
            fontsize=16
        )
    
        ax.set_xlim(0, 20)                     # existing line
        ax.xaxis.set_major_locator(mtick.MaxNLocator(integer=True))
        ax.tick_params(axis='x', labelsize=16)
    
    
        # draw bin edges
        
        edges = sorted({b[0] for b in bins} | {b[1] + 1 for b in bins})
        for edge in range(0, 251, 5):      # 0,5,10,15,…,250
            if edge==0:
                ax.axvline(edge, linestyle="-", color='black', linewidth=2,label='bins')
            else:
                ax.axvline(edge, linestyle="-", color='black', linewidth=2)
                
    
        # mark observation
        ax.axvline(x=obs, color='red', linestyle='--',label='Number of reported cases', linewidth=5)
        # ax.text(0.81, 0.95, f"log score = {score:.2f}",
        #         transform=ax.transAxes, verticalalignment="top", fontsize=16)
    
        # optional legend (comment out if undesired)
        ax.legend(loc='upper right', fontsize=15)
        plt.savefig(f'Logscore_{excluded_year}_{day/7}_orange-1year.jpg', format='jpg', bbox_inches='tight')
        plt.savefig(f'Logscore_{excluded_year}_{day/7}_orange-1year.pdf', format='pdf', bbox_inches='tight')
    
        plt.show()
    #print("scores",log_scores)
    # ---- summary printout -----------------------------------------------
    
    
    
    # for (day, obs, lam), sc in zip(MSE_list, log_scores):
    #     print(f"Day {day}: obs={obs}, λ={lam:.2f}, log score={sc:.2f}")
    
    
    
    
    
    
    
    
    
    #print("mean of historical data",mean_r_excluding(excluded_year, folder_path="."))
    
    
    
    
    rate_fixed=(mean_r_excluding(excluded_year, folder_path="."))
    
    
    
    
    
    log_scores_null = []
    d=[]
    lam=rate_fixed
    
    
    for week_num, (day, obs, lam) in enumerate(MSE_list, 1):
        idx = bin_index_for_value(obs)
        if idx is None:
            score = -10.0
        else:
            lam=rate_fixed
            P = bin_probability(lam, idx)
            score = np.log(P) if P > 0 else -10.0
            score = max(score, -10.0)
        log_scores_null.append(score)
        d.append(int(day/7))
        #print(d)
    
        # ----------------- plotting --------------------------------------
        k_vals  = np.arange(0, 20)
        pmf_vals = poisson_pmf(k_vals, lam)
    
        # fig, ax = plt.subplots(figsize=(14, 8))
        # ax.plot(k_vals, pmf_vals, '--', linewidth=4)
    
        # # --- NEW: shade the bin containing the observation ---------------
        # if idx is not None:
        #     lo, hi = bins[idx]
        #     mask = (k_vals >= lo) & (k_vals <= hi)
        #     ax.fill_between(k_vals[mask], 0, pmf_vals[mask],
        #                     color='green', alpha=0.35, label=f'-log(Observed‐bin area)={score:.2f}', )
    
        # ax.set_xlabel("Number of cases", fontsize=16)
        # ax.set_ylabel("$P_{Null}(k | λ)$", fontsize=16)
        # ax.set_title(
        #     f"Week {day/7} ({excluded_year}): Null Poisson Distribution (λ ={rate_fixed}) "
        #     "for Orange County, California",
        #     fontsize=16
        # )
    
        # ax.set_xlim(0, 20)                     # existing line
        # ax.xaxis.set_major_locator(mtick.MaxNLocator(integer=True))
        # ax.tick_params(axis='x', labelsize=16)
    
    
        # # draw bin edges
        
        # edges = sorted({b[0] for b in bins} | {b[1] + 1 for b in bins})
        # for edge in range(0, 251, 5):      # 0,5,10,15,…,250
        #     if edge==0:
        #         ax.axvline(edge, linestyle="-", color='black', linewidth=2,label='bins')
        #     else:
        #         ax.axvline(edge, linestyle="-", color='black', linewidth=2)
                
    
        # # mark observation
        # ax.axvline(x=obs, color='red', linestyle='--', label='Number of reported cases',linewidth=5)
        # # ax.text(0.81, 0.95, f"log score = {score:.2f}",
        # #         transform=ax.transAxes, verticalalignment="top", fontsize=16)
    
        # # optional legend (comment out if undesired)
        # ax.legend(loc='upper right', fontsize=15)
    
        # plt.show()
        
        
        
                
        
    T_log_scores.append([sum(log_scores),sum(log_scores_null)])    
    count = sum(1 for a, b in zip(log_scores, log_scores_null) if a > b)
    em=count+em
    tm=len(log_scores_null)+tm
    #print("number ", count)  # also outputs: 2
    
        
        
        
        
        
        
        
        

        
        
        
        
        
        
        
        
        
        
    
    # ---- summary printout -----------------------------------------------
    # for (day, obs, lam), sc in zip(MSE_list, log_scores_null):
        
        
    #     print(f"Day {day}: obs={obs}, λ={lam:.2f}, log score={sc:.2f}")
    
    
    
    
    #print("DDDD",d)
    

    
    plt.figure(figsize=(14, 8))
    plt.plot(d, log_scores, 'k--*', label='Predicted Poisson process', linewidth=5, markersize=14)
    plt.plot(d, log_scores_null, 'g-*', label='Null model', linewidth=5, markersize=14)

    plt.xlabel('Number of weeks with reported cases', fontsize=16)
    plt.ylabel('Logarithmic Score', fontsize=16)
    plt.title(f'Comparison Between Null Model and Poisson Model for Spillover Severity in Orange County, California (Excluding {excluded_year})', fontsize=16)
    plt.legend(fontsize=15)
    
    plt.xticks(d, fontsize=14)  # ✅ Only show actual x values from d
    plt.yticks(fontsize=14)

    plt.tight_layout()
    plt.savefig(f'comparison_plot_{excluded_year}_orange-1year.jpg', format='jpg', bbox_inches='tight')
    plt.show()
    











arr = np.array(T_log_scores)
arrmse = np.array(comparison_MSE)

#print("log",arr)

# First bar plot: Log Scores
fig, ax = plt.subplots(figsize=(10, 5))
bar_width = 0.35
x = np.arange(len(Y))

ax.bar(x - bar_width/2, arr[:, 0], bar_width, label='Model Sum(log_scores)')
ax.bar(x + bar_width/2, arr[:, 1], bar_width, label='Null Sum(log_scores_null)')
ax.set_xlabel('Year')
ax.set_ylabel('Annual sum of Log_Score')
ax.set_title('Time Series of Total Log Scores for Each Year Prediction - Orange County')
ax.set_xticks(x)
ax.set_xticklabels(Y, rotation=45)
ax.legend()
plt.tight_layout()
plt.savefig('comparison_plot_Total_log_score_orange_bar.jpg', format='jpg', bbox_inches='tight')
plt.show()

# Second bar plot: RMSE
fig, ax = plt.subplots(figsize=(10, 5))
ax.bar(x - bar_width/2, arrmse[:, 0], bar_width, label='Model NRMSE')
ax.bar(x + bar_width/2, arrmse[:, 1], bar_width, label='Null NRMSE')
ax.set_xlabel('Year')
ax.set_ylabel('Annual NRMSE')
ax.set_title('Time Series of NRMSEs for Each Year Prediction - Orange County')
ax.set_xticks(x)
ax.set_xticklabels(Y, rotation=45)
ax.legend()
plt.tight_layout()
plt.savefig('comparison_plot_Total_RMSE_orange_bar.jpg', format='jpg', bbox_inches='tight')
plt.show()











# print("**********************",em)
# print(point_estimation)





# Calculate absolute errors
errors = [abs(obs - pred) for obs, pred in point_estimation]

# Total number of predictions
total = len(errors)
# print("number of all predictions",total)

# Cumulative thresholds
thresholds = [1, 2, 3, 4, 5]

# Initialize cumulative counts
cumulative_counts = {f"<{t}": 0 for t in thresholds}
cumulative_counts[">=5"] = 0  # changed '>5' to '>=5' for clarity

# Count how many errors fall under each threshold
for e in errors:
    for t in thresholds:
        if e < t:
            cumulative_counts[f"<{t}"] += 1
    if e >= 5:
        cumulative_counts[">=5"] += 1

# Convert to percentage
percentages = {k: round((v / total) * 100, 2) for k, v in cumulative_counts.items()}
# print(percentages)



# print("Rmse",arrmse[:, 0])

# #print("Rmse",plot_excluded_year_kde_curve_and_daily_risk(YEARS, excluded_year, C_predict_file='C_predict.txt')[1])

# print("MSE",Comaoprison_realMSE)





