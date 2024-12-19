# PredSim gait conditions

This repository contains the source code for the research article [ref]. In this study we compared measured and simulated kinematics, kinetics, ground reaction forces and metabolic power in various gait conditions.

### dependencies

This code depends on several software packages

1. And old (and slightly adapted) version of PredSim ([GitHub - KULeuvenNeuromechanics/PredSim: Generate predictive simulations of human locomotion](https://github.com/KULeuvenNeuromechanics/PredSim)) is that basis of this repository. This old version is already included in the repository (folder PredSim). Please note that for future simulations projects I advice you to work with the most recent version of PredSim (and not this one).

2. Casadi is used to formulate the NLP ( https://web.casadi.org/)

3. Matlab opensim api (https://opensimconfluence.atlassian.net/wiki/spaces/OpenSim/pages/53089380/Scripting+with+Matlab)

4. Visual studio (to compile opensimAD code)

### installation instructions

Installation

1. clone/fork this repository

2. add repository to your matlab path

3. Install compiler code to create .dll files from osim models. Do this in matlab by running the script: ./PredSim/InstallOsim2Dll. This software is based on [GitHub - antoinefalisse/opensimAD: Libraries for OpenSimAD - OpenSim with support for Algorithmic Differentiation.](https://github.com/antoinefalisse/opensimAD)

4. Add Neuromechanics toolkit as a submodule *git submodule init*

## replicate simulations

**1. convert .osim models for specific studies**

go to the folder ./PredSim/AdaptOpenSimModel to convert your opensim model for a specific study. For example the script (ConvertModelsBrowning.m) adds mass to the opensim model in the location specified in Browning et al. 2008. This scripts saves the .osim file, .cpp and .cdll files to solve inverse dynamics (needed for implicit formulation of skeleton dynamics) in the folders Subjects. You probably have to adapt the main path to the repository and the visual studio compiler (e.g. Visual Studio 17 2022)

**2. Run predictive simulation**

go to the folder ./PredSim and simulate a specific experiment. For example the script Simulate_Browning2008 predicts the walking motion as in the experiment of Browning.

### Single-muscle simulations

The single muscle simulations to evaluate the mechanical efficiency of the muscle+metabolic energy model can be found in *PredSim\Efficiency_muscle_model*. For example the script *....\MaxEfficiencySimModel.m* generates figure 4 in the paper. The script *....\MaxEfficiencySimModel_varEModels.m* generates figure 14 in the appendix.

### Compare to experimental data

You can compare the measured and simulated step frequency and metabolic cost using the script *.../PlotFigures/TableSimResults.m.* This script generates figures 1,2,3 and 5. This scripts assumes that you replicated all experiments in simulations. If you did this, you will probably have to adapt some directories to specify where you saved the simulations results. If you want to use my simulations results you can download them using this link and point to this folder in the script TableSimResults.

I you want a detailed comparison with the experiment of McDonald 2021 you can run the script *.../PlotFigures/DetailedFigureMcDonald.m*

I you want to detailed comparison with the dataset of van der Zee 2022 with various walking speeds you can run the script *.../PlotFigures/DetailedAnalysis_EffectWalkingSpeed.m*. Note that the experimental data from vanderZee was processed with Addbiomechanics (https://addbiomechanics.org/download_data.html) and gait cycle average data was computed with a python script (*AddBiomechanics\src\Kinematics_Kinetics_VanDerZee.py*) using the addbiomechanics API ([Working with AddBiomechanics Data &mdash; Nimble Physics 0.4.0 documentation](https://nimblephysics.org/docs/working-with-addbiomechanics-data.html)).

### download simulation results and processed experimental data

If you don't want to convert all osim models to .dll files, rerun all simulations and process data with addbiomechanics (but just want to play with the sim results) you can download data from zenodo: https://doi.org/10.5281/zenodo.14524619

Make sure that you point in your matlab code to these folder where needed (and add the subjects folder to ./PredSim/Subjects) if you to recreate my figures using my code.
