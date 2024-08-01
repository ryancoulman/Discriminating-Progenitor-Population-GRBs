import pandas as pd
import matplotlib.pyplot as plt
import numpy as np
from scipy.optimize import curve_fit
import matplotlib as mpl
from astropy.modeling.models import BrokenPowerLaw1D
from astropy.io import fits
import requests
from bs4 import BeautifulSoup
import re
import os
from joblib import Parallel, delayed
from typing import Optional, Dict, Any



#/////////////// Functions to read in lc files, and GRB durations /////////////////////



# just converts all the names of the grb lc files in a given folder into an array to later run through
def GRBS_in_folder(folder_name, file_suffix, grbname_or_fname='grb'):

    # allows sorting by date 
    def extract_date(s):
        return s[3:8]

    fname_list = [file for file in os.listdir(folder_name) if file.endswith(file_suffix)]
    fname_list = sorted(fname_list, key=extract_date)

    grbname_list = []
    for file_name in fname_list:
        grb_name = file_name.split('_')[0]
        grbname_list.append(grb_name)

    return grbname_list, fname_list


# Reads in fermi lc files into pandas dataframe 
def read_fermi_fits_file(filename, Tmin, Tmax, const=1000):

    # provide full path to file 
    folder = filename.split('/')[0]
    grbname = (filename.split('/')[1]).split('_')[0]
    # i had a problem manipulating the fermi files read in initially (something about a big endian error) and so this creates a csv file of the fermi data 
    #   instead of  afits file. if created before then will read that in
    if os.path.exists(f'{folder}/{grbname}_copy.csv'):
        new_df = pd.read_csv(f'{folder}/{grbname}_copy.csv')
    else:
        lcdata = fits.getdata(f'{filename}', 1)
        df = pd.DataFrame(lcdata)

        new_df = pd.DataFrame(columns=['TIME', 'RATE'])
        i = 0
        for t, r in zip(df['Time'], df['Rate']):
            new_df.loc[i] = [t, r]
            i += 1

        new_df.to_csv(f'{folder}/{grbname}_copy.csv', index=False)

    # drops all data outside tmin and tmax 
    new_df = new_df[(new_df['TIME'] > Tmin) & (new_df['TIME'] < Tmax)]
    new_df.reset_index(drop=True, inplace=True)

    # add constant to all data (cant log negatives) and log 
    const = max(abs(new_df['RATE'].min())+10, const)
    new_df['log_15-350'] = np.log(new_df['RATE'] + const)

    return new_df

# returns tmin and tmax from provided dataframe provided column names align 
def get_Tmin_Tmax(grb, df):

    if df['GRB'].isin([grb]).any():
        GRB_row = df[df['GRB'] == grb]
    else:
        print(f"The GRB {grb} is not in the DataFrame {df}")

    Tmin = GRB_row['Tmin'].iloc[0]
    Tmax = GRB_row['Tmax'].iloc[0]
    if 'Terr' in df.columns:
        Terr = GRB_row['Terr'].iloc[0]
    else:
        Terr = pd.NA

    return Tmin, Tmax, Terr

# trivial function to print the grb name currently being evaluated 
def print_banner(name, character="/", length=60):
    banner = character * length + " " + name + " " + character * length
    print(banner)


# used to plot either the graph of Sigma or the LC depending upon the kwarg
def plot_func(data, x_axis, y_axis, Sig_or_LC, window_size=1, tpre=None, t90=None, yerror=None):
    mpl.rcParams['font.size'] = 15

    if Sig_or_LC == 'Sig':

        plt.xscale('log')
        plt.yscale('log')
        plt.xlabel('Δt (s)')
        plt.ylabel('$\sigma_{\Delta t}$')

    elif Sig_or_LC == 'LC':

        if window_size > 1:
            GRB_average = data.rolling(window=window_size).mean().iloc[window_size - 1::window_size].reset_index(
                drop=True)
            data = GRB_average
        if tpre is not None:
            for t in tpre:
                plt.axvline(x=t, color='r', linestyle='--')
        if t90 is not None:
            for t in t90:
                plt.axvline(x=t, color='g', linestyle='--')

        plt.xlabel('Time since BAT / seconds')
        plt.ylabel('Count rate')

    if yerror is None:
        plt.plot(data[f'{x_axis}'], data[f'{y_axis}'], marker='o', linestyle='', color='black', markersize=2)
    else:
        plt.errorbar(data[f'{x_axis}'], data[f'{y_axis}'], yerr=data[f'{yerror}'], marker='o', linestyle='', color='black', markersize=2, capsize=2)
    

    return plt.show()



# ///////////////// Functions to calcuate the allan variance and MVT //////////////////////////



# calcuate the allan varaiance (sigma_deltat)
# uses joblib package to compute on parralel cores. turn fast to false if causing problems (only results in small speed boost unless propogating errors too)
def calc_sigma(curve, bin_width, scale, fast=True):

    def process_scale(curve, bin_width, scale):
        x_bar_col = curve['log_15-350'].rolling(window=scale, min_periods=1).mean().iloc[scale - 1::scale].reset_index(drop=True)
        haar_coefficients_sq_col = (x_bar_col.diff(periods=1).iloc[1::2].reset_index(drop=True).values) ** 2
        sum_of_haar_sq = np.sum(haar_coefficients_sq_col)
        sigma = np.sqrt(sum_of_haar_sq * (scale / len(curve)))
        return {'delta_t': scale * bin_width / 1000.0, 'sigma': sigma}

    scale_range = int(len(curve) / (2 * scale)) 
    bin_scales = [i * scale for i in range(1, scale_range + 1)]

    if fast:
        results = Parallel(n_jobs=-1)(delayed(process_scale)(curve, bin_width, scale) for scale in bin_scales)
    else:
        results = [process_scale(curve, bin_width, scale) for scale in bin_scales]

    data = pd.DataFrame(results)

    return data


# average the sigma values in enlarging intervals to maintain uniform spacing in log space 
def actual_avg(sigma_values, initial_interval=1, multiplier=np.e):
    
    num_points = len(sigma_values)
    
    # Initialize variables
    interval = initial_interval
    start_index = 0
    averaged_values = []
    
    # Loop until the end of data
    while start_index < num_points:
        # Calculate end index for current interval
        end_index = min(start_index + interval, num_points)
        
        # Calculate average for current interval
        average = sigma_values.iloc[int(start_index):int(end_index)].mean()
        averaged_values.append(average)
        
        # Update start index for next interval
        start_index = end_index
        
        # Update interval size
        interval = interval * multiplier
        
        # Ensure interval is at least 1
        interval = max(interval, 1)
    
    # Create a DataFrame with averaged values
    averaged_df = pd.DataFrame(averaged_values, columns=sigma_values.columns)
    averaged_df.columns = ['delta_t', 'sigma_avg']
    
    return averaged_df


# determines the roots of the polynomial 
def classify_root(roots, poly, min_time, max_time):

    local_minima = []
    local_maxima = []

    roots = roots[(roots > min_time) & (roots < max_time)]

    for root in roots:

        value = 10 ** poly(np.log10(root))
        prev_value = 10 ** poly(np.log10(root - 1e-5))
        next_value = 10 ** poly(np.log10(root + 1e-5))

        if prev_value > value and next_value > value:
            local_minima.append(root)
        elif prev_value < value and next_value < value:
            local_maxima.append(root)

    if len(local_minima) > 1:
        local_minima = [min(local_minima)]
    if len(local_maxima) > 1:
        local_maxima = [max(local_maxima)]

    return local_minima, local_maxima


# fits a polynomial to the data and determines roots using classify_roots() to find the local minima and maxima 
def fit_poly(data, x_axis, y_axis, degree=7):

    x_axis = np.array(data[f'{x_axis}'])
    y_axis = np.array(data[f'{y_axis}'])

    min_time = min(x_axis)
    max_time = max(x_axis)

    # Fit a polynomail to data
    coefficients = np.polyfit(np.log10(x_axis), np.log10(y_axis), degree)
    poly = np.poly1d(coefficients)

    # Find roots in log scale then convert back to linear
    derivative_poly = poly.deriv()
    roots_log = derivative_poly.roots
    real_roots = 10 ** roots_log.real

    local_minima, local_maxima = classify_root(real_roots, poly, min_time, max_time)


    if len(local_minima) >= 1 and len(local_maxima) >= 1:
        local_minima = local_minima[0]
        local_maxima = local_maxima[0]
    elif len(local_minima) == 0 and len(local_maxima) > 0:
        local_minima = min(x_axis)
        local_maxima = local_maxima[0]
        print("No min found")
    elif len(local_maxima) == 0 and len(local_minima) > 0:
        local_minima = local_minima[0]
        local_maxima = max(x_axis)
        print("No max found")
    else:
        local_minima = min(x_axis)
        local_maxima = max(x_axis)

    if local_minima > local_maxima:
        local_maxima = max(x_axis)

    return local_minima, local_maxima, poly


# Plot the second derivative of the polynomial and find the roots corresponding to the maximum of the first derivative polynomial (the MVT) 
# Also determines the error as the differnece in MVT between the degree above and below 
def plot_second_derivative_errors(data, x_axis, y_axis, local_minima, local_maxima, poly, degree, yerror=None):

    x_axis = np.array(data[f'{x_axis}'])
    y_axis = np.array(data[f'{y_axis}'])
    if yerror is not None:
        yerror = np.array(data[f'{yerror}'])

         
    coefficients_higher = np.polyfit(np.log10(x_axis), np.log10(y_axis), degree+1)
    poly_deghigher = np.poly1d(coefficients_higher)
    coefficients_lower = np.polyfit(np.log10(x_axis), np.log10(y_axis), degree-1)
    poly_deglower = np.poly1d(coefficients_lower)

    poly_array = [poly_deglower, poly_deghigher, poly]
    # Find roots in log scale then convert back to linear

    mvt = []

    for i, pol in enumerate(poly_array):

        derivative_poly = pol.deriv()
        second_derivative_poly = derivative_poly.deriv()

        # Find max of first derivative between local min and max
        roots_log = second_derivative_poly.roots
        real_roots = 10 ** roots_log.real
        try:
            filtered_roots = real_roots[(real_roots > local_minima) & (real_roots < local_maxima)]
            filtered_roots = sorted(filtered_roots)

            length = len(filtered_roots)
            if length != 0:
                first_root_val = 1
                i = 0
                while first_root_val > 1e-6 and i <= length - 1:
                    first_root_val = np.abs(second_derivative_poly(np.log10(filtered_roots[i])))
                    first_deriv_max = filtered_roots[i]
                    i += 1
            else:
                print(f"No roots in range under 1e-6 for {'deglow' if i == 0 else 'deghigh' if i == 1 else 'poly deg'}")
                first_deriv_max = np.nan

        except Exception:
            filtered_roots = np.nan
            print(f"Could not find any roots in range for {'deglow' if i == 0 else 'deghigh' if i == 1 else 'poly deg'}")
            first_deriv_max = np.nan

        mvt.append(first_deriv_max)

    mvt_error = ( (abs(mvt[2]-mvt[0]) if not np.isnan(mvt[0]) else 0) + (abs(mvt[2]-mvt[1]) if not np.isnan(mvt[1]) else 0) ) / 2
    mvt_error = 0 if np.isnan(mvt_error) else mvt_error
    first_deriv_max = mvt[2]

    # Plot the polynomial, first derivative, and second derivative
    mpl.rcParams['font.size'] = 15
    plt.figure(figsize=(8, 6))
    # Plot the polynomial
    if yerror is not None: 
        plt.errorbar(x_axis, y_axis, yerr=yerror, marker='o', linestyle='', color='black', markersize=2, capsize=2)
    else:
        plt.plot(x_axis, y_axis, marker='o', linestyle='', color='black', markersize=2)
    x_values_poly = np.linspace(min(x_axis), max(x_axis), 50000)
    plt.plot(x_values_poly, 10 ** poly(np.log10(x_values_poly)), color='red')
    plt.scatter(first_deriv_max, 10 ** poly(np.log10(first_deriv_max)), color='purple', label='first deriv max')
    plt.scatter(first_deriv_max + mvt_error, 10 ** poly(np.log10(first_deriv_max + mvt_error)), color='blue', label='mvt err')
    plt.scatter(first_deriv_max - mvt_error, 10 ** poly(np.log10(first_deriv_max  - mvt_error)), color='blue', label='mvt err')
    plt.scatter(local_maxima, 10 ** poly(np.log10(local_maxima)), color='green', label='Local Maxima')
    plt.scatter(local_minima, 10 ** poly(np.log10(local_minima)), color='orange', label='Local Minima')
    plt.xscale('log')
    plt.yscale('log')
    plt.xlabel('$\Delta t$ (s)')
    plt.ylabel('$\sigma_{Δt}$')
    plt.legend()
    plt.show()
 
    plt.figure(figsize=(12, 8))
    # Plot the first derivative
    plt.subplot(3, 1, 1)
    plt.plot(x_values_poly, derivative_poly(np.log10(x_values_poly)), color='green', label='First Derivative')
    plt.scatter(first_deriv_max, derivative_poly(np.log10(first_deriv_max)), color='purple', label='first deriv max')
    plt.scatter(local_maxima, derivative_poly(np.log10(local_maxima)), color='green')
    plt.scatter(local_minima, derivative_poly(np.log10(local_minima)), color='orange')
    plt.axhline(0, color='black', linestyle='--', linewidth=1)
    plt.xscale('log')
    plt.ylabel('First Derivative')
    plt.legend()

    # Plot the second derivative
    plt.subplot(3, 1, 2)
    y_values_second_derivative = second_derivative_poly(np.log10(x_values_poly))
    plt.plot(x_values_poly, y_values_second_derivative, color='purple', label='Second Derivative')
    plt.scatter(filtered_roots, second_derivative_poly(np.log10(filtered_roots)), color='cyan',
                label='other roots in range')
    plt.scatter(first_deriv_max, second_derivative_poly(np.log10(first_deriv_max)), color='purple',
                label='first deriv max')
    plt.scatter(local_maxima, second_derivative_poly(np.log10(local_maxima)), color='green')
    plt.scatter(local_minima, second_derivative_poly(np.log10(local_minima)), color='orange')
    plt.axhline(0, color='black', linestyle='--', linewidth=1)
    plt.xscale('log')
    plt.xlabel('Δt / seconds')
    plt.ylabel('Second Derivative')
    plt.legend()

    plt.tight_layout()
    plt.show()

    print(f"mvt: {first_deriv_max} +- {mvt_error}")
    return first_deriv_max, mvt_error

# //////////////// classes to help with analysis ///////////////



# dataframe to store mvt products during analysis 
class MVTDataFrame:
    def __init__(self, filename=None):
        self.filename = filename 
        self.columns = {'GRB': str, 'T90_duration': float, 'T90_duration_err': float, 'T90_MVT': float, 'T90_MVT_err': float, 'Tpre_MVT': float, 'Tpre_MVT_err': float, 'Tpre_duration': float}
        if filename and os.path.exists(f'MVTs/{filename}'):
            self.retrieve_df()
        else:
            self.df = pd.DataFrame(columns=self.columns.keys()).astype(self.columns)
        self.col_keys = list(self.columns.keys())

    def add_row(self, row, grb=None, t90_duration=None, t90_dur_err=None, tpre_duration=None, t90_mvt=None, tpre_mvt=None, col_name=None, t90_mvt_err=None, tpre_mvt_err=None):
            if grb is not None:
                self.df.loc[row, 'GRB'] = grb
            if t90_duration is not None:
                self.df.loc[row, 'T90_duration'] = t90_duration
            if tpre_duration is not None:
                self.df.loc[row, 'Tpre_duration'] = tpre_duration
            if t90_mvt is not None:
                self.df.loc[row, 'T90_MVT'] = t90_mvt
            if tpre_mvt is not None:
                self.df.loc[row, 'Tpre_MVT'] = tpre_mvt
            if col_name is not None:
                if col_name[0] == 'T90':
                    self.df.loc[row, 'T90_MVT'] = col_name[1]
                    self.df.loc[row, 'T90_MVT_err'] = col_name[2]
                    self.df.loc[row, 'T90_duration'] = col_name[3]
                    self.df.loc[row, 'T90_duration_err'] = col_name[4] if len(col_name) > 4 else pd.NA
                elif col_name[0] == 'Tpre':
                    self.df.loc[row, 'Tpre_MVT'] = col_name[1]
                    self.df.loc[row, 'Tpre_MVT_err'] = col_name[2]
                    self.df.loc[row, 'Tpre_duration'] = col_name[3]
            if t90_mvt_err is not None:
                self.df.loc[row, 'T90_MVT_err'] = t90_mvt_err
            if tpre_mvt_err is not None:
                self.df.loc[row, 'Tpre_MVT_err'] = tpre_mvt_err
            if t90_dur_err is not None:
                self.df.loc[row, 'T90_duration_err'] = t90_dur_err
            self.save_df()

    def get_dataframe(self):
        return self.df

    def get_column_keys(self):
        return self.col_keys

    def save_df(self):
        self.df.to_csv(f'MVTs/{self.filename}', index=False)

    def retrieve_df(self):
        self.df = pd.read_csv(f'MVTs/{self.filename}')
        self.df = self.df.astype(self.columns)

    def get_grbs_done(self):
        return self.df['GRB'].values 
    


# not all grbs provide a clear allan varainace and mvt. This checks if you are happy with the mvt automatically calcuated and allows you to re-run if not
# note this does not work well in the native jupyter notebook due to the nature of input popups. isntead run in vscode for a better experience 
class Refine:

    def __init__(self, Tmin:float, Tmax:float,  scale:int, degree:int, multiplier:float):
        self.scale: int = scale
        self.degree: int = degree
        self.local_min:  Optional[float] = None
        self.local_max:  Optional[float] = None
        self.Tmin: float = Tmin
        self.Tmax: float = Tmax
        self.multiplier: float = multiplier
        self.reset = 'n'
        self.initial_values = self._get_initial_values()

    def repeat(self, row, mvtClass, Tx='T90'):
            answer = input("Does this look reasonable (y/n/nan/c_{val}): ")
            if answer == 'y':
                return False
            elif answer == 'nan':
                if Tx == 'T90':
                    mvtClass.add_row(row=row, t90_mvt=pd.NA, t90_mvt_err=pd.NA)
                elif Tx == 'Tpre':
                    mvtClass.add_row(row=row, tpre_mvt=pd.NA, tpre_mvt_err=pd.NA)
                return False
            elif 'c_' in answer:
                value = float(answer.split('_')[1])
                if Tx == 'T90':
                    mvtClass.add_row(row=row, t90_mvt=value)
                elif Tx == 'Tpre':
                    mvtClass.add_row(row=row, tpre_mvt=value)
                return False
            else:
                self.update_defaults()
                return True

    def update_defaults(self):
        for key, value in self.__dict__.items():
            param = input(f"Enter new value for {key} (press Enter to keep current value): ")
            if param != '':
                if key in ['scale', 'degree']:  # Specify keys that should be integers
                    setattr(self, key, int(param))
                elif key in ['local_min', 'local_max', 'Tmin', 'Tmax', 'multiplier']:  # Specify keys that should be floats
                    setattr(self, key, float(param))
                else:
                    setattr(self, key, param)
            if self.reset == 'y':
                self.reset_defaults()

    def assign_stationairy(self, local_minima, local_maxima):
        if self.local_min is None:
            self.local_min = local_minima
        if self.local_max is None:
            self.local_max = local_maxima

    def _get_initial_values(self):
        return {
            'scale': self.scale,
            'degree': self.degree,
            'multiplier': self.multiplier,
            'local_min': self.local_min,
            'local_max': self.local_max,
            'Tmin': self.Tmin,
            'Tmax': self.Tmax,
            'multiplier': self.multiplier,
            'reset': self.reset
        }

    def reset_defaults(self):
        for key, value in self.initial_values.items():
            setattr(self, key, value)

