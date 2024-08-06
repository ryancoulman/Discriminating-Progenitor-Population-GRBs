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

# *Project Abstract*

*The recent detection of the long-duration, merger originating GRB211211A and GRB230307A has sparked a renewed interest in means to discriminate the progenitor population of GRBs based on their gamma-ray light curve alone. Traditional classification schemes, reliant solely on temporal and spectral properties to group GRBs into long and short bursts associated with stellar collapse and compact binary mergers, respectively, have proven inadequate in accurately discerning these anomalous GRBs. Building upon previous work, we utilise Haar wavelets to measure the variability of the bulk emission and systematically apply the technique in a novel manner to pulses preceding the main event, termed precursors. We find the former analysis of the prompt emission unable to disseminate the anomalous population, whereas the latter study of precursors correctly demonstrates clustering of the anomalies with the short merger population. Therefore, by utilising precursors, our work aids in the accurate classification of GRBs, which is essential for leveraging these bursts to understand star formation, metal enrichment, and the universe’s evolution.*

## Features

- **Object-Orientatied Programming :** Calculates Bessel functions of the first kind $J_{m}(x)$ for orders $m = 0$ to $m = 9$.
- **Error handling:** Generates high-precision values suitable for scientific and engineering applications.
- **Tabular Data Output:** Outputs results in a tabular format for easy analysis and interpretation.

## Key Concepts Demonstrated

- **Object-Orientatied Programming:** Utilises functional C programming for mathematical calculations.
- **Error Handling:** Implements the CTR to approximate the integral to a high degree of accuracy using $N = 10,000$ sub-intervals.
- **Handling Data-Sets:** using **Pandas** 
- **Pandas, Numpy, Matplotlib, Scipy:** Develops efficient algorithms to perform complex mathematical computations.
- **Web-Scraping:** Once calcuated, samples the Bessel function at varied m and x values using a smaller N to verify the precision to the number of decimal places rquested by the user, improving computational efficiency by reducing N when possible if the output is unchanged.
- **Parralel-Core Processing:**


## Code Overview

### Main Program (`main.c`)

- Defines and initializes necessary variables and constants.
- Computes Bessel functions for orders $m = 0$ to $m = 9$.
- Outputs the results in a structured tabular format to `Bessel_Output.txt`.

### Output File (`Bessel_Output.txt`)

- Contains computed Bessel function values for $x$ ranging from 0 to 99.
- Organized in a table where each row corresponds to a different $x$ value, and each column corresponds to a Bessel function of a particular order.

## Results

inlcude conclusin in report where say what achived 
atttache most important plot of mvt90 vs mvt pre
