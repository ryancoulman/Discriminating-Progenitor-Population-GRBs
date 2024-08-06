# Overview

The aim of my final year university project was to computationally classify the progenitor population of long and short gamma-ray bursts (GRBs) using their lightcurve (LC) alone. This involved developing a pipeline that leverages data analysis libraries in Python to process and analyse GRB data.

I used `pandas` dataframes to store the LC flux data in order to compute a mathmatical function termed a Haar scalorgarm,

$$
\sigma_{\Delta t}^{2} = \frac{\Delta t}{t}\sum\limits_{i=0}^{t/2\Delta t -1}h^{2}_{i, \Delta t},
$$

given,

$$
h_{i, \Delta t} = \overline{X}_{2i+1, \Delta t} - \overline{X} _{2i, \Delta t},
$$

and where $\overline{X}$ is the average natural logarithm of the flux over ∆t consecutive bins, $X$. From this I used `numpy` to fit an $n^{th}$ order polynomial to the data, which typically follows a sigmoidal function, and determined its first derivative maximum between the local minima and maxima to define the minimum variable timescale (MVT) of interest. Short of identifying the best fit polynomial visually, I set up a fully automatic pipeline to compute the MVT of around 100 GRBs.  

The early work analysing the MVT presented in this project provides a promising avenue for discriminating the progenitor source of GRBs, and is currently undergoing further validation for potential publication. This project also received a mark of 80%. 

## Key Concepts Demonstrated

- **Data Manipulation with Pandas:** Utilises pandas for efficient handling, filtering, and transformation of large datasets (up to 1 million rows).
- **UNIX Command-Line Tools:** Frequent use of UNIX commands to download GRBs and implemented a Bash script to automate the process for a array of GRB names. 
- **Scientific Computing with NumPy:** Employs NumPy for numerical calculations, array operations, and fitting of polynomials.
- **Matplotlib for Data Visualization:**  Uses Matplotlib extensively for plotting data.
- **Parallel Processing:** Implements joblib's Parallel library for parallel computations, significantly reducing processing time.
- **File Handling and OS Operations:** Deploys the `os` module for handling file operations and navigating the file system.
- **Web Scraping with BeautifulSoup:** Uses requests and BeautifulSoup for extracting data from the *Swift* database to automatuically extract the duration of a given GRB 
- **Polynomial Fitting and Derivatives:** Fits polynomial models to data and analyzes their derivatives to find critical points.
- **Error Propagation:** Calculates and includes error propagation in Haar Scalogram calculations.
- **Data Export:** Manages data export to csv, txt, and xlsx files for later use or analysis.
- **Statistical Analysis:** Conducts statistical analysis of findings using the spearmanr class from `SciPy`
- **Astropy:** Uses Astropy for handling FITS files and performing specific astronomical calculations.

## Code Overview

### Neccassary set of functions and classes (`mvt_functions.py`)

- Contains all functions essential for reading light curve files, calculating the MVT, and processing data
- Reading and sorting GRB light curve files.
- Extracting and processing GRB names and dates.
- Error handling and parallel processing using joblib.
- Implementing statistical models and fitting data using scipy and astropy.

## Main program (`main.ipynb`)

- A Jupyter notebook that calls functions from mvt_functions.py to process GRB data.
- Data import and cleaning.
- Matching and merging data from different GRB databases.
- Calculations of various time-related parameters for GRBs.
- Plotting and visualising MVT results using matplotlib.

## Data analaysis of the MVT (plot_data.ipynb`)

- Another Jupyter notebook focused on processing the results obrtained from main.ipynb
- Detailed data processing and merging steps from independent GRB samples.
- Calculation of statistical correlations and fits.
- Advanced plotting and visualisation techniques.
  
## Results

inlcude conclusin in report where say what achived 
atttache most important plot of mvt90 vs mvt pre

--- 

# *Project Abstract*

*The recent detection of the long-duration, merger originating GRB211211A and GRB230307A has sparked a renewed interest in means to discriminate the progenitor population of GRBs based on their gamma-ray light curve alone. Traditional classification schemes, reliant solely on temporal and spectral properties to group GRBs into long and short bursts associated with stellar collapse and compact binary mergers, respectively, have proven inadequate in accurately discerning these anomalous GRBs. Building upon previous work, we utilise Haar wavelets to measure the variability of the bulk emission and systematically apply the technique in a novel manner to pulses preceding the main event, termed precursors. We find the former analysis of the prompt emission unable to disseminate the anomalous population, whereas the latter study of precursors correctly demonstrates clustering of the anomalies with the short merger population. Therefore, by utilising precursors, our work aids in the accurate classification of GRBs, which is essential for leveraging these bursts to understand star formation, metal enrichment, and the universe’s evolution.*
