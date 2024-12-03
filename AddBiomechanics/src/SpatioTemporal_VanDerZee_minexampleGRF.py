import nimblephysics as nimble
from typing import List
import numpy as np
import matplotlib.pyplot as plt
import time
from utilities import *

# minimal example to show problem with zero GRF at first frame

datapath = "/home/maarten/Documents/data/addbiomech/VanDerZee2022"
# Trials to analyze (we want free walking at various speeds
trial_names_sel = np.array([2])
nsubj = 1
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
        grfMat_raw = np.zeros([nfr_tot, 6])
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
                        grf = frame_to_render.processingPasses[iproc].groundContactForce
                        grfMat[ifr, :] = grf
                        # get raw forces ()
                        rawGRF = frame_to_render.rawForcePlateForces
                        grfMat_raw[ifr, 0:3] = rawGRF[0]
                        grfMat_raw[ifr, 3:6] = rawGRF[1]
                        tMat[ifr] = ifr * your_subject.getTrialTimestep(trial)
                        ifr = ifr + 1
                        if frame == 0:
                            print(grf) # filtered grf is always zero at the first frame
                            print(rawGRF) # raw grf not :)


        # compute gait events for all trials
        tMat = np.delete(tMat, np.s_[ifr:nfr_tot], 0)
        grfMat = np.delete(grfMat, np.s_[ifr:nfr_tot], 0)
        grfMat_raw = np.delete(grfMat_raw, np.s_[ifr:nfr_tot], 0)
        [ths_l, tto_l, hs_l, to_l] = detect_heelstrike_toeoff(tMat, grfMat[:, 1], 100, dtOffPlate=0.2)
        [ths_r, tto_r, hs_r, to_r] = detect_heelstrike_toeoff(tMat, grfMat[:, 1], 100, dtOffPlate=0.2)
        stridefreq  = 1/((np.nanmedian(np.diff(ths_l)) + np.nanmedian(np.diff(ths_r)))/2)

        # norm stride freq to subject size
        stridefreq_norm[s,i] = stridefreq/(your_subject.getHeightM()*0.5)

        plt.figure()
        plt.subplot(2,1,1)
        plt.plot(grfMat[:,1])
        plt.plot(grfMat[:,4])
        plt.xlabel('Frames')
        plt.ylabel('Force-filt [N]')
        plt.subplot(2,1,2)
        plt.plot(grfMat_raw[:,1])
        plt.plot(grfMat_raw[:,4])
        plt.xlabel('Frames')
        plt.ylabel('Force [N]')
        plt.legend(['Fy1','Fy2'])

plt.show()




