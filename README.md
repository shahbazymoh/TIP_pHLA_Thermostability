**A pipeline for thermostability profiling of HLA-bound peptides** 

**Thermal immunopeptidome profiler (TIP) code, version 1.2** 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%% This code has been created, developed, and validated by Mohammad (Moh) Shahbazy 

%%% A PhD candidate at the Department of Biochemistry and Molecular Biology 

%%% Clayton Campus, Monash University (Melbourne, Australia)

%%% Laboratory: Purcell Lab, Biomedicine Discovery Institute (BDI) 

%%% Acknowledgements: Dr Nathan Croft and Prof Anthony Purcell, for supervision

%%% Copyright© Mohammad Shahbazy (2019-present).


%%% Creation date: June 2019

%%% Last modification: March 2023

% EXAMPLE DATASET: C1R-B*57:01 data acquired by Fusion Orbitrap (Thermo)

% DIA-NN version 1.8.0 was used for library-based DIA data processing

% Compatible with the next versions of DIA-NN (e.g., 1.8.1 or 1.9.2)

%%%%%%%%%%%%%%%*************** IMPORTANT ***************%%%%%%%%%%%%%%%%%%%

% Input files: Please search your DIA data by DIA-NN (*MBR option should be unchecked) 
% and use the "report.tsv" to organize and manipulate the results in a new "csv" file 
% according to the example "data_table" below (imported data in the 1st step)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% This TIP code (Thermal immunopeptidome profiler) performs thermostability profiling for individual HLA-bound peptides. 
% The following steps are included in the code to preprocess and manipulate the DIA-based quantitative immunopeptidomics 
dataset to calculate melting temperature (Tm) for peptide-HLA complexes (pHLAs).

%%% *Step 1) Importing Data

%%% *Step 2) Importing the target list of the peptide sequences

%%% *Step 3) selecting the quantification level for profiling

%%% *Step 4) Manipulating data to preprocess and generate a "data matrix" of the quantified precursors across data points

%%% *Step 5) Isolating iRT peptides' data for the normalization to iRTs

%%% *Step 6) Rolling the data points to smooth the profile and manipulating outliers with constant values

%%% *Step 7) Removing precursors with missing values at any first three datapoints of this dataset (37C to 46C)

%%% *Step 8) Sigmoid curve fitting to estimate Tm value (IC50) for thermostability profiling of immunopeptidomes

%%% *Step 9) Filtering valid Tm values based on the expected range and removing failed fittings 

%%% *Step 10) Evaluation/Validation of the fitting quality 
% Percentage of the fitted curves with a good correlation with experimental data points, R >= 0.75  

%%% *Step 11) Basic visualization to check results: Plotting a histogram to show the distribution for the calculated Tm values

%%% *Step 12) Drawing Denaturation Profiles with the corresponding Tm values for individual HLA peptides

%%% *** Export the results for downstream immunoinformatics and other analyses
******************************************************************************************************

**Contact** For technical support or any questions, please send an email to the following developers:

Dr Mohammad Shahbazy: shahbazymoh@gmail.com or mohammad.shahbazy@monash.edu

Dr Nathan Croft: nathan.croft@monash.edu
