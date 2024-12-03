
import nimblephysics as nimble
from typing import List
import numpy as np
import matplotlib.pyplot as plt
import time
from utilities import ButterFilter_Low_NaNs
# Test with addbiomechanics data
#--------------------------------

datapath = "/home/maarten/Documents/data/addbiomech/tutorial"
geometrypath = "/home/maarten/Documents/data/opensim/geometry/"
# Load the model
your_subject = nimble.biomechanics.SubjectOnDisk(datapath + "/Falisse2017_subject_1.b3d")
# Read the skeleton that was optimized by the first process pass (always kinematics)
# Use the geometryFolder argument to specify where to load the bone mesh geometry from
skeleton: nimble.dynamics.Skeleton = your_subject.readSkel(
    processingPass=0,
    geometryFolder=geometrypath)


# A B3D file contains multiple trials of data, so we need to select which one we want to load
ntrials = your_subject.getNumTrials()

for trial in range(0,ntrials):

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
    ndof = your_subject.getNumDofs()
    qstore = np.zeros([nfr, ndof])
    tau_store = np.zeros([nfr, ndof])
    grf_store = np.zeros([nfr, 9])
    time = np.zeros(nfr)
    ekin_store = np.zeros(nfr)

    for frame in range(0, nfr):
        # Get the frame we want to render
        frame_to_render = trial_frames[frame]

        # Set the skeleton's state to the state in the frame
        skeleton.setPositions(frame_to_render.processingPasses[iproc].pos)
        skeleton.setVelocities(frame_to_render.processingPasses[iproc].vel)

        # get external forces
        cop = frame_to_render.processingPasses[iproc].groundContactCenterOfPressure
        grf = frame_to_render.processingPasses[iproc].groundContactForce
        grf_t = frame_to_render.processingPasses[iproc].groundContactTorque

        # get the inverse dynamics (I assume that his is a "top down" approach without external forces)
        skeleton.getInverseDynamics(frame_to_render.processingPasses[iproc].acc)

        # store position information
        qstore[frame, :] = frame_to_render.processingPasses[iproc].pos
        tau_store[frame, :] = frame_to_render.processingPasses[iproc].tau # is the joint moment ?
        time[frame] = frame_to_render.t * your_subject.getTrialTimestep(trial)
        grf_store[frame, :] = grf
        ekin_store[frame] = skeleton.computeKineticEnergy()

    # test lowpass filter on processed data (this is not ideal but should give a very similar result)
    fs = np.round(1./your_subject.getTrialTimestep(trial))
    tau_filtered = ButterFilter_Low_NaNs(fs, tau_store)


    # try to plot a figure -- kinematics generalized coords
    # plt.figure()
    # dofs = skeleton.getDofs()
    # for n in range(0, ndof):
    #     ax = plt.subplot(4, 6, n+1)
    #     ax.plot(time, qstore[:, n])
    #     ax.set_title(dofs[n].getName())
    #
    # try to plot a figure-- moments
    plt.figure()
    dofs = skeleton.getDofs()
    for n in range(0, ndof):
        ax = plt.subplot(4, 6, n+1)
        ax.plot(time, tau_store[:, n])
        ax.plot(time, tau_filtered[:, n])
        ax.set_title(dofs[n].getName())

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


