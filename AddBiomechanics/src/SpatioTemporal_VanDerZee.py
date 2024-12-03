import nimblephysics as nimble
from typing import List
import numpy as np
import matplotlib.pyplot as plt
import time
from utilities import *


datapath = "/home/maarten/Documents/data/addbiomech/VanDerZee2022"
# Trials to analyze (we want free walking at various speeds
trial_names_sel = np.array([2, 5, 8, 11, 14, 16, 32])
trail_speeds = np.array([0.7, 0.9, 1.1, 1.6, 1.8, 2, 1.4])
nsubj = 10
# pre allocate output
stridefreq_norm = np.zeros([nsubj+1, len(trial_names_sel)])

# loop over all subjects
for s in range(1, nsubj+1):

    subjectfolder = "p" + str(s)
    filename = "p" +str(s) + ".b3d"

    # Load the model
    filepath = datapath + "/" + subjectfolder + "/" + filename
    your_subject = nimble.biomechanics.SubjectOnDisk(filepath)

    # Load all the frames from the selected trial
    nTrials = your_subject.getNumTrials()
    iproc = 0 # raw data
    # find all indices of the selected trials
    nfr_tot = 2000*6
    for i in range(0, len(trial_names_sel)):
        #get the selected trial
        iTrial = trial_names_sel[i]
        iTrialStart = 'trial' + str(iTrial) + '_s'
        # pre allocate
        grfMat = np.zeros([nfr_tot, 6])
        tMat = np.zeros([nfr_tot])
        ifr = 0
        # loop over all trials
        for trial in range(0, nTrials):
            tname = your_subject.getTrialName(trial)
            tname_sel = tname[0:len(iTrialStart)]
            if (iTrialStart == tname_sel):
                # read all the frames
                nfr = your_subject.getTrialLength(trial)
                trial_frames: List[nimble.biomechanics.Frame] = your_subject.readFrames(
                    trial=trial,
                    startFrame=0,
                    numFramesToRead=your_subject.getTrialLength(trial))
                if nfr > 100:
                    for frame in range(0, nfr):
                        frame_to_render = trial_frames[frame]
                        #grf = frame_to_render.processingPasses[iproc].groundContactForce
                        #grfMat[ifr, :] = grf
                        rawGRF = frame_to_render.rawForcePlateForces
                        grfMat[ifr, 0:3] = rawGRF[0]
                        grfMat[ifr, 3:6] = rawGRF[1]
                        tMat[ifr] = ifr * your_subject.getTrialTimestep(trial)
                        ifr = ifr + 1

        # compute gait events for all trials
        tMat = np.delete(tMat, np.s_[ifr:nfr_tot], 0)
        grfMat = np.delete(grfMat, np.s_[ifr:nfr_tot], 0)
        [ths_l, tto_l, hs_l, to_l] = detect_heelstrike_toeoff(tMat, grfMat[:, 1], 100, dtOffPlate=0.2)
        [ths_r, tto_r, hs_r, to_r] = detect_heelstrike_toeoff(tMat, grfMat[:, 4], 100, dtOffPlate=0.2)
        stridefreq  = 1/((np.nanmedian(np.diff(ths_l)) + np.nanmedian(np.diff(ths_r)))/2)

        # norm stride freq to subject size
        stridefreq_norm[s,i] = stridefreq/(your_subject.getHeightM()*0.5)



    print('finished with subject ' + str(s))

print(stridefreq_norm)
# save the stride frequency file to a .csv file
fileout = datapath + '/output/' + 'stridefreq.csv'
with open(fileout, "w", newline="") as f:
    writer = csv.writer(f)
    # write the header
    header = ["subject", "walkspeed", "stridefreq"]
    writer.writerow(header)
    # write the data
    for s in range(1, nsubj+1):
        for i in range(0, len(trial_names_sel)):
            data = [s, trail_speeds[i], stridefreq_norm[s,i]]
            writer.writerow(data)

print('finished with program')







