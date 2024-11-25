# PredSim gait conditions

This repository contains the source code for the research article [ref]. In this study we ....

### dependencies

This code depends on several software packages

1. And old (and slightly adapted) version of PredSim ([GitHub - KULeuvenNeuromechanics/PredSim: Generate predictive simulations of human locomotion](https://github.com/KULeuvenNeuromechanics/PredSim)) is that basis of this repository. This old version is already included in the repository (folder PredSim). 

2. Casadi is used to formulate the NLP (ref casadi)

3. Matlab opensim api (ref)

### installation instructions

Installation

1. clone/fork this repository

2. Install compiler code to create .dll files from osim models. Do this in matlab by running the script: ./PredSim/InstallOsim2Dll.



## replicate simulations

**1. convert .osim models for specific studies**

go to the folder ./PredSim/AdaptOpenSimModel to convert your opensim model for a specific study. For example the script (ConvertModelsBrowning.m) adds mass to the opensim model in the location specified in Browning et al. (ref).

**2. Run predictive simulation**

go to the folder ./PredSim and simulate a specific experiment. For example the script Simulate_Browning2008 predicts the walking motion as in the experiment of Browning.



### Single-muscle simulations

The single muscle simulations to evaluate the mechanical efficiency of the muscle+metabolic energy model can be found in *PredSim\Efficiency_muscle_model*. For example the script *....\MaxEfficiencySimModel.m* generates figure 4 in the paper. The script *....\MaxEfficiencySimModel_varEModels.m* generates the figure 14 in the appendix.

### Compare to experimental data

You can compare the measured and simulated step frequency and metabolic cost using the script *./PlotFigures/TableSimResults.m.*. This script generates figures 1,2,3 and 5. 




