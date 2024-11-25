ReadMe
_________________________________________________________________________________________________________________________
Author:		Kirsty A. McDonald
Date: 		August, 2022
Software:	Mathworks MATLAB
_____________________________________________________________________________________________________________________
Processed dataset used for: 
		McDonald K.A., Cusumano, J.P., Hieronymi, A., Rubenson, J.R.
		Humans trade-off whole-body energy cost to avoid overburdening muscles while walking 
____________________________________________________________________________________________________________________

FILE: EnergyCost_FatigueAvoidance_MasterData.mat
 
Variables within .mat file are described below:
____________________________________________________________________________________________________________________

SubjectParams (structure):

	
	level one	PWS		preferred walking speed (m/s)
			Age		years
			Gender		f=female; m=male
			BodyMass	body mass (kg)
			Height		height (m)
			Transition	% incline transitioned from choosing incline to choosing crouching (first incline value that crouch was chosen as the preferred locomotor state)
			[data: ten subjects, rows (if 10x1) or columns (if 1x10) represent subjects (S01-S10; subject ID)]
____________________________________________________________________________________________________________________

C_metP, C_asq, C_avol, C_amax:

	10x6 array	rows		participants (N=10)
			cols		mean crouch, 0% incline, 6% incline, 12% incline, 18% incline, 24% incline

Individual processing methods for these variables are described in the corresponding manuscript.

	