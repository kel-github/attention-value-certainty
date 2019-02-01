# Incentive value and spatial certainty combine additively to determine visual priorities

This repos contains the stimulation code and analysis code/data for the 2 experiments presented in Garner, Bowman & Raymond, Incentive value and spatial certainty combine addditively to determine visual priorities
BioRxiv: 


(c) Kelly Garner, 2019

Any questions contact:
getkellygarner@gmail.com

### Acknowledgment

You are allowed to use this data and software for free, but please acknowledge/cite the authors if you publish your efforts


### how to use this repository:
For both experiments the following structure applies:
~~~~
		code-analysis-and-task/
					/Contains EXP1/ and EXP2/ (See below)
					/Exp_1_2_ANALYSIS.R is the code for the main analysis (presented in the results section), with some additional NHST analyses. Notes for running are documented at the top of the code 
					/Exp_1_2_Figures.R produces figures 2 and 3 from the paper. Again, see notes at the top of the code for the guide to running 
					/EXP[1 or 2]
				   			/TASK  Contains task/stimulation code. See notes at the top of the run code - e.g. cue_probs_v*.m for dependencies and instructions for running
							EXP2 has extra .mat files - these contain variables for the predetermined counterbalancing across participants
					/wideform_summary_data/ contains group summary data files in wideform, for RT and Accuracy respectively, so that non-R users can recreate the analyses in the paper in their own software packages
							legend:
							.4/.6 | .2/.8 etc = cue certainty condition
							valid | invalid   = cue validity
							high | low        = value of the location where the target appeared
							static | decay (exp 2 only) = reward condition
~~~~