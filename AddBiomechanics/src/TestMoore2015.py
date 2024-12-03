
import nimblephysics as nimble
from typing import List
import numpy as np
import matplotlib.pyplot as plt
import time
from utilities import *
# Test with addbiomechanics data
#--------------------------------

# one of the main problems here is that the trail tags are not used. I'm missing the treadmill speed here and the
# perturbed and unperturbed trials

datapath = "/home/maarten/Documents/data/addbiomech/Moore2015"
geometrypath = "/home/maarten/Documents/data/opensim/geometry/"
subjectfolder = "subject12"
filename = "subject12.b3d"

# Load the model
your_subject = nimble.biomechanics.SubjectOnDisk(datapath + "/" + subjectfolder + "/" + filename)

# Read the skeleton that was optimized by the first process pass (always kinematics)
# Use the geometryFolder argument to specify where to load the bone mesh geometry from
skeleton: nimble.dynamics.Skeleton = your_subject.readSkel(
    processingPass=0,
    geometryFolder=geometrypath)


# A B3D file contains multiple trials of data, so we need to select which one we want to load
#ntrials = your_subject.getNumTrials()


# pre allocate the output matrices
nfr_tot = 2000*40
ndof = 23
qMat = np.zeros([nfr_tot, ndof])
qdMat = np.zeros([nfr_tot, ndof])
tauMat = np.zeros([nfr_tot, ndof])
tMat = np.zeros([nfr_tot])
grfMat = np.zeros([nfr_tot, 6])
ifr = 0



for trial in range(1, 33):

   # Load all the frames from the selected trial
    trial_frames: List[nimble.biomechanics.Frame] = your_subject.readFrames(
        trial=trial,
        startFrame=0,
        numFramesToRead=your_subject.getTrialLength(trial))

    # get GRF bodies
    grf_bodies = your_subject.getGroundForceBodies()

    # Loop through all the frames, and render them to the GUI
    #frame = 0
    iproc = 2
    nfr = your_subject.getTrialLength(trial)
    if nfr > 100:
        ndof = your_subject.getNumDofs()
        #qstore = np.zeros([nfr, ndof])
        #qdotstore = np.zeros([nfr, ndof])
        #tau_store = np.zeros([nfr, ndof])
        #grf_store = np.zeros([nfr, 6])
        #time = np.zeros(nfr)
        #ekin_store = np.zeros(nfr)

        for frame in range(0, nfr):
            # Get the frame we want to render
            frame_to_render = trial_frames[frame]

            # Set the skeleton's state to the state in the frame
            if len(frame_to_render.processingPasses)>iproc:
                skeleton.setPositions(frame_to_render.processingPasses[iproc].pos)
                skeleton.setVelocities(frame_to_render.processingPasses[iproc].vel)

                # get external forces
                cop = frame_to_render.processingPasses[iproc].groundContactCenterOfPressure
                grf = frame_to_render.processingPasses[iproc].groundContactForce
                grf_t = frame_to_render.processingPasses[iproc].groundContactTorque

                # get the inverse dynamics (I assume that his is a "top down" approach without external forces)
                skeleton.getInverseDynamics(frame_to_render.processingPasses[iproc].acc)

                # store position information
                qMat[ifr, :] = frame_to_render.processingPasses[iproc].pos
                qdMat[ifr, :] = frame_to_render.processingPasses[iproc].vel
                tauMat[ifr, :] = frame_to_render.processingPasses[iproc].tau # is the joint moment ?
                tMat[ifr] = ifr* your_subject.getTrialTimestep(trial)
                grfMat[ifr, :] = grf
                #ekin_store[ifr] = skeleton.computeKineticEnergy()
                ifr = ifr+1

# remove not used pre-allocated rows
qMat= np.delete(qMat, np.s_[ifr:nfr_tot], 0)
qdMat= np.delete(qdMat, np.s_[ifr:nfr_tot], 0)
tauMat= np.delete(tauMat, np.s_[ifr:nfr_tot], 0)
tMat= np.delete(tMat, np.s_[ifr:nfr_tot], 0)
grfMat = np.delete(grfMat,np.s_[ifr:nfr_tot], 0)

# test lowpass filter on processed data (this is not ideal but should give a very similar result)
fs = np.round(1./your_subject.getTrialTimestep(trial))
tau_filtered = ButterFilter_Low_NaNs(fs, tauMat)
qdot_filtered = ButterFilter_Low_NaNs(fs, qdMat)
JointPower = tau_filtered * qdot_filtered

# test some event detection on this data
[ths, tto, hs, to]= detect_heelstrike_toeoff(tMat, grfMat[:,1], 50, dtOffPlate=0.2)

# plot a more generic figure
plt.figure()
ct = 0
idofs_plot = [6, 9, 10]
dofs = skeleton.getDofs()
for n in idofs_plot:
    ax = plt.subplot(4, 3, ct + 1)
    ax.plot(tMat, qMat[:, n])
    ax.set_title(dofs[n].getName())

    ax = plt.subplot(4, 3, ct + 1 + 3)
    ax.plot(tMat, tau_filtered[:, n])

    ax = plt.subplot(4, 3, ct + 1 + 6)
    ax.plot(tMat, grfMat[:, ct])
    if ct == 1:
        plt.vlines(ths,0,800,colors='k', linestyles='solid')

    ax = plt.subplot(4, 3, ct + 1 + 9)
    ax.plot(tMat, JointPower[:, n])
    ct = ct + 1

# test plot gait cycle average data
# dt_stride_min = np.nanmean(np.diff(ths))
[grfMat_cycle, DatMean, DatSTD, DatMedian] = norm_cycle(tMat, ths, grfMat, treshold_dtStride_max=1.9,
                                                 treshold_dtStride_min=1, bool_print = True,
                                                 n_points_int = 100)

[qMat_cycle, DatMean, DatSTD, DatMedian] = norm_cycle(tMat, ths, qMat, treshold_dtStride_max=1.9,
                                                 treshold_dtStride_min=1, bool_print = True,
                                                 n_points_int = 100)

plt.figure()
for ind in range(0, 3):
    ax = plt.subplot(2, 3, ind+1)
    plt.plot(grfMat_cycle[:, ind,:])
    ax = plt.subplot(2, 3, ind+4)
    plt.plot(qMat_cycle[:, idofs_plot[ind],:])

    # try to plot a figure -- kinematics generalized coords
    #plt.figure()
    #dofs = skeleton.getDofs()
    #for n in range(0, ndof):
    #    ax = plt.subplot(4, 6, n+1)
    #    ax.plot(time, qstore[:, n])
    #    ax.set_title(dofs[n].getName())

    # try to plot a figure-- moments
    # plt.figure()
    # dofs = skeleton.getDofs()
    # for n in range(0, ndof):
    #     ax = plt.subplot(4, 6, n+1)
    #     ax.plot(time, tau_store[:, n])
    #     ax.plot(time, tau_filtered[:, n])
    #     ax.set_title(dofs[n].getName())

    # try to plot a figure -- GRF
    # plt.figure()
    # for n in range(0, 9):
    #     ax = plt.subplot(3, 3, n+1)
    #     ax.plot(time, grf_store[:, n])
    #
    # plt.figure()
    # plt.plot(time, ekin_store)
    # plt.show()
plt.show()


