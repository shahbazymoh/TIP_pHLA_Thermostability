%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
**A pipeline for thermostability profiling of HLA-bound peptides** 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
**Themal immunopeptidome profiler (TIP) code, version 1.2** 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%% This code has been created, developed, and validated by Mohammad (Moh) Shahbazy 
%%% A PhD canidate at Department of Biochemistry and Molecular Biology  
%%% Clayton Campus, Monash University (Melbourne, Australia)
%%% Laboratory: Purcell Lab, Biomedicine Discovery Institute (BDI) 
%%% Acknowledgements: Dr Nathan Croft and Prof Anthony Purcell, for supervison
%%% Copyright© Mohammad Shahbazy (2019-2023).


%%% Creation date: June 2019
%%% Last modification: March 2023

% EXAMPLE DATASET: C1R-B*57:01 data acquired by Fusion Orbitrap (Thermo)

% DIA-NN version 1.8.0 was used for library-based DIA data processing
% Compatible with the next versions of DIA-NN (e.g., 1.8.1 or 1.9.2)

%%%%%%%%%%%%%%%*************** IMPORTANT ***************%%%%%%%%%%%%%%%%%%%
% Input files: Please search your DIA data by DIA-NN (*MBR option should be uncheck) 
% and use the "report.tsv" to organize and manupulate the results in a new "csv" file 
% according to the example "data_table" below (imported data in the 1st step)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


% This TIP code (Themal immunopeptidome profiler) performs thermostability profiling for individual HLA-bound peptides. 
% The following steps are included in the code to preprocess and manupulate the DIA-based quantitative immunopeptidomics 
dataset to calculate melting temperature (Tm) for peptide-HLA complexes (pHLAs).


%%% *Step 1) Importing Data

%%% *Step 2) Importing target list of the peptide sequences

%%% *Step 3) selecting the quantification level for profiling

%%% *Step 4) Manupulating data to preprocess and generate a "data matrix" of the quantified precursors across data points

%%% *Step 5) Isolating iRT peptides' data for the normalization to iRTs

%%% *Step 6) Rolling the data points to smooth the profile and manupulating outliers with constant values

%%% *Step 7) Removing precursors with missing values at the any first three datapoints of this dataset (37C to 46C)

%%% *Step 8) Sigmoid curve fitting to estimate Tm value (IC50) for thermostability profiling of immunopeptidomes

%%% *Step 9) Filtering valid Tm values based on the expected range and removing failed fittings 

%%% *Step 10) Evaluation/Validation of the fitting quality 
% Percentage of the fitted curves with a good correlation with experimental data points R >= 0.75  

%%% *Step 11) Basic visualization to check results: Plotting a histogram to show the distribution for the calculated Tm values

%%% *Step 12) Drawing Denaturation Profiles with the correspong Tm values for individual HLA peptides


******************************************************************************************************

Contact For technical support or any questions, please send an email to the following developers:
Dr Mohammad Shahbazy: shahbazymoh@gmail.com or mohammad.shahbazy@monash.edu
Dr Nathan Croft: nathan.croft@monash.edu
