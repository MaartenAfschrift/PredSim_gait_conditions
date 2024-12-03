# Create .csv file with kinematics and kinetics of ankle, knee and hip angle

import nimblephysics as nimble
from typing import List
import numpy as np
import matplotlib.pyplot as plt
import time
from utilities import *
import pandas as pd


geometrypath = "/home/maarten/Documents/data/opensim/geometry/"
datapath = "/home/maarten/Documents/data/addbiomech/VanDerZee2022"

# Trials to analyze (we want free walking at various speeds
trial_names_sel = np.array([2, 5, 8, 11, 14, 16, 32])
trail_speeds = np.array([0.7, 0.9, 1.1, 1.6, 1.8, 2, 1.4])
nsubj = 10
# pre allocate output
stridefreq_norm = np.zeros([nsubj+1, len(trial_names_sel)])

# loop over all subjects
IK_LargeMat = np.full([100, 23, nsubj+1, len(trial_names_sel)], np.nan)
ID_LargeMat = np.full([100, 23, nsubj+1, len(trial_names_sel)], np.nan)
GRF_LargeMat = np.full([100, 6, nsubj+1, len(trial_names_sel)], np.nan)
for s in range(1, nsubj+1):

    subjectfolder = "p" + str(s)
    filename = "p" +str(s) + ".b3d"

    # Load the model
    filepath = datapath + "/" + subjectfolder + "/" + filename
    your_subject = nimble.biomechanics.SubjectOnDisk(filepath)

    skeleton: nimble.dynamics.Skeleton = your_subject.readSkel(
        processingPass=0,
        geometryFolder=geometrypath)

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
        ndof = your_subject.getNumDofs()
        IKMat = np.zeros([nfr_tot,ndof])
        IDMat = np.zeros([nfr_tot,ndof])
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
                        # get ground reaction forces
                        rawGRF = frame_to_render.rawForcePlateForces
                        grfMat_raw[ifr, 0:3] = rawGRF[0]
                        grfMat_raw[ifr, 3:6] = rawGRF[1]
                        grfMat[ifr,:] = frame_to_render.processingPasses[iproc].groundContactForce
                        # update time vector
                        tMat[ifr] = ifr * your_subject.getTrialTimestep(trial)
                        # get joint kinematics
                        # Set the skeleton's state to the state in the frame
                        #skeleton.setPositions(frame_to_render.processingPasses[iproc].pos)
                        #skeleton.setVelocities(frame_to_render.processingPasses[iproc].vel)
                        IKMat[ifr, :] = frame_to_render.processingPasses[2].pos
                        IDMat[ifr, :] = frame_to_render.processingPasses[2].tau  # is the joint moment ?
                        ifr = ifr + 1

        # compute gait events for all trials
        tMat = np.delete(tMat, np.s_[ifr:nfr_tot], 0)
        grfMat = np.delete(grfMat, np.s_[ifr:nfr_tot], 0)
        grfMat_raw = np.delete(grfMat_raw, np.s_[ifr:nfr_tot], 0)
        IKMat = np.delete(IKMat, np.s_[ifr:nfr_tot], 0)
        IDMat = np.delete(IDMat, np.s_[ifr:nfr_tot], 0)
        [ths_l, tto_l, hs_l, to_l] = detect_heelstrike_toeoff(tMat, grfMat[:, 1], 100, dtOffPlate=0.2)
        [ths_r, tto_r, hs_r, to_r] = detect_heelstrike_toeoff(tMat, grfMat[:, 4], 100, dtOffPlate=0.2)

        # problem with ground reaction forces
        iallzeros = np.where((grfMat[:, 1] == 0) & (grfMat[:, 4] == 0))
        grfMat[iallzeros, : ] = np.nan

        # filter kinematics and kinetics
        fs = np.round(1. / your_subject.getTrialTimestep(0))
        IDMat = ButterFilter_Low_NaNs(fs, IDMat)
        IKMat = ButterFilter_Low_NaNs(fs, IKMat)
        IDMat = IDMat/ your_subject.getMassKg()
        grfMat = grfMat/ your_subject.getMassKg()
        #grfMat = ButterFilter_Low_NaNs(fs, grfMat)

        # norm kinematics and kinetics to stride
        mean_cycle = np.nanmean(np.diff(ths_l))
        [IK_all, IK_mean, IK_std, IK_median] = norm_cycle(tMat, ths_l, IKMat,
                                                          treshold_dtStride_max = mean_cycle + 0.2,
                                                          treshold_dtStride_min=mean_cycle - 0.2)
        [ID_all, ID_mean, ID_std, IDmedian] = norm_cycle(tMat, ths_l, IDMat,
                                                         treshold_dtStride_max=mean_cycle + 0.2,
                                                         treshold_dtStride_min=mean_cycle - 0.2)
        [GRF_all, GRF_mean, GRF_std, GRFmedian] = norm_cycle(tMat, ths_l, grfMat,
                                                             treshold_dtStride_max=mean_cycle + 0.2,
                                                             treshold_dtStride_min=mean_cycle - 0.2)

        # we want to print the median IK and ID data to a .csv file
        headers = []
        dofs  = skeleton.getDofs()
        for ih in range(0,your_subject.getNumDofs()):
            headers.append(dofs[ih].getName())

        speed_str = str(int(np.round(trail_speeds[i]*10)))
        filename_IK = datapath + '/output/subject_' + str(s) + '_' + speed_str + '_IK.csv'
        filename_ID = datapath + '/output/subject_' + str(s) + '_' + speed_str + '_ID.csv'
        filename_GRF = datapath + '/output/subject_' + str(s) + '_' + speed_str + '_GRF.csv'

        df = pd.DataFrame(IK_median, columns=headers)
        df.to_csv(filename_IK)

        df = pd.DataFrame(IDmedian, columns=headers)
        df.to_csv(filename_ID)

        grf_headers = ['Flx','Fly','Flz','Frx','Fry','Frz']
        df = pd.DataFrame(GRFmedian, columns=grf_headers)
        df.to_csv(filename_GRF)

        # store in large matrix
        IK_LargeMat[:, :, s-1, i] = IK_median
        ID_LargeMat[:, :, s-1, i] = IDmedian
        GRF_LargeMat[:, :, s-1, i] = GRFmedian

    print('finished with subject ' + str(s))

 # compute mean for each walking speed
IK_mean = np.nanmean(IK_LargeMat,axis = 2)
ID_mean = np.nanmean(ID_LargeMat,axis = 2)
GRF_mean = np.nanmean(GRF_LargeMat,axis = 2)
IK_std = np.nanstd(IK_LargeMat,axis = 2)
ID_std = np.nanstd(ID_LargeMat,axis = 2)
GRF_std = np.nanstd(GRF_LargeMat,axis = 2)
for i in range(0, len(trial_names_sel)):
    speed_str = str(int(np.round(trail_speeds[i] * 10)))
    filename_IK = datapath + '/output/mean_' + speed_str + '_IK.csv'
    filename_ID = datapath + '/output/mean_' + speed_str + '_ID.csv'
    filename_GRF = datapath + '/output/mean_' + speed_str + '_GRF.csv'
    df = pd.DataFrame(np.squeeze(IK_mean[:, :, i]), columns=headers)
    df.to_csv(filename_IK)
    df = pd.DataFrame(np.squeeze(ID_mean[:, : , i]), columns=headers)
    df.to_csv(filename_ID)
    df = pd.DataFrame(np.squeeze(GRF_mean[:, :, i]), columns=grf_headers)
    df.to_csv(filename_GRF)

    filename_IK = datapath + '/output/mean_' + speed_str + '_IK_std.csv'
    filename_ID = datapath + '/output/mean_' + speed_str + '_ID_std.csv'
    filename_GRF = datapath + '/output/mean_' + speed_str + '_GRF_std.csv'
    df = pd.DataFrame(np.squeeze(IK_std[:, :, i]), columns=headers)
    df.to_csv(filename_IK)
    df = pd.DataFrame(np.squeeze(ID_std[:, : , i]), columns=headers)
    df.to_csv(filename_ID)
    df = pd.DataFrame(np.squeeze(GRF_std[:, :, i]), columns=grf_headers)
    df.to_csv(filename_GRF)



        # plt.figure()
        # dofs = skeleton.getDofs()
        # for n in range(0, ndof):
        #     ax = plt.subplot(4, 6, n + 1)
        #     ax.plot(IK_all[:, n, :],'b')
        #     ax.plot(IK_mean[:, n],'k')
        #     ax.set_title(dofs[n].getName())
        #
        # plt.figure()
        # for n in range(0, 6):
        #     ax = plt.subplot(2, 3, n + 1)
        #     ax.plot(GRF_all[:, n, :],'b')
        #     ax.plot(GRF_mean[:, n],'k')
        #
        # #del IKMat, IDMat, grfMat
        #
        # plt.figure()
        # plt.plot(tMat, IKMat[:,13])
        # plt.vlines(ths_l,-1, 1,'k')
        #
        # plt.figure()
        # plt.plot(tMat, IKMat[:,6])
        # plt.vlines(ths_r,-1, 1,'k')
        # #
        # plt.figure()
        # plt.plot(tMat, grfMat[:,1])
        # plt.vlines(ths_l,-1, 500,'k')
        # #
        # plt.figure()
        # plt.plot(tMat, grfMat[:,1],'b')
        # plt.plot(tMat, grfMat_raw[:, 4], 'r')



