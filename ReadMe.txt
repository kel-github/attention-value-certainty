REPOSITORY NOTES:

This repos contains the stimulation code and analysis code/data for the 2 experiments presented in Garner & Raymond, Distinct selection mechanisms for when predictions and rewards guide visual selective attention
https://docs.google.com/document/d/1j_9irYF9LbegaxuxnWxSFdmkLe8nSrm6AsWl-ShVCYI/edit?usp=sharing

Any questions contact:
getkellygarner@gmail.com

For both experiments the following structure applies
ADDBIAS_REPOS/
Contains EXP1/ and EXP2/ (See below)
Exp_1_2_ANALYSIS.R is the code for the main analysis (presented in the results section), with some additional NHST analyses - conducted out of curiosity to see if it corroborated the Bayes analysis.
Exp_1_2_Figures.R produces figures 2 and 3 from the paper.

EXP[1 or 2]/

TASK/ 
Contains task/stimulation code
% DEPENDENCIES
% coded on Matlab 2013a using PTB v3.0.12
% run on Stone SOFREP-144
% with Asus VG278HE
EXP2 has extra .mat files - these contain variables for the predetermined counterbalancing across participants

ANALYSIS/
.csv files = raw data
The remaining code is a mixture of .R, .Rmd, and .html files - I have been playing with formats for communicating results hence the mixture.
The steps outlined in the ‘data cleaning’ section of the paper are carried out in the .Rmd files. The .html files reflect the code and the plotted data (e.g. so you can see the piecewise regression fits). The output of this is the .R data file called something along the lines of exp#_clean_BS_fv_date.R
Any bespoke functions called by the .Rmd file will be in the folder as an .R file with lower case lettering in the title
The file named EXP#_ANALYSIS_MIXDMDLS_BIAS.R is the workspace file, which contains all the objects created during the analysis (see Exp_1_2_ANALYSIS.R)

Confusing things - sometimes files refer to exps 3 and 5, instead of 1 and 2 - this is because of versions I piloted that I have not added to the repository.