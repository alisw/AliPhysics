// $Id$

/*! 

\page README_calign Offline Calibration and Alignment
 
\section calign_s1 Offline monitoring of the Alignment and Calibration

The macro MUONClusterInfo.C reads the ESDs or RecPoints and creates a tree of 
AliMUONClusterInfo objects which can be conveniently used to plot track to 
cluster residuals, cluster and pads properties.

Usage:
<pre>
root [0] .L $ALICE_ROOT/MUON/MUONClusterInfo.C+
root [1] MUONClusterInfo(nevents, esdFileName, inFileName, outFileName); 
</pre>

 1) <code> if ( esdFileName != "" ) </code>: 

 Loop over ESD event and fill AliMUONClusterInfo object with cluster + 
 corresponding track parameters.

 2) <code> if ( inFileName != "" ) </code>:
  
 Loop over RecPoints and fill AliMUONClusterInfo object with cluster not 
 attached to a track;

 3) Write results in a new root file.

\section calign_s2 How to run offline calibration and alignment tasks

AliMUONAlignmentTask:
 
It is an AliAnalysisTask to align the MUON spectrometer. The Task reads 
as input ESDs and feeds the MUONTracks to AliMUONAlignment. The alignment 
itself is performed by AliMillepede. A OCDB entry is written with the 
alignment parameters. A list of graph are written to monitor the alignment 
parameters. 

AliMUONReAlignTask:

It is an AliAnalysisTask to realign the MUON spectrometer. The Task reads 
as input ESDs moves the clusters of a MUONTrack acoording to the re aligned 
geometry taken from a misalignment file in the OCDB and refit the track. 
Then it produces a AliMUONClusterInfo object for each cluster. The output is 
a TTree of AliMUONClusterInfo.

AliMUONChamberCalibrationTask:

It is a AliAnalysisTaskSE to recalibrate the MUON spectrometers. The output is 
a TTree of AliMUONClusterInfo. 


The AliMUON*Task are designed to run in an analyis train of the Analysis 
Framework. The macros AddTaskMuon*.C can be used to add the task to an 
analysis train. An usage example based on the new standard AnalysisTrainNew.C 
will be provided soon. 


This chapter is defined in the READMEcalign.txt file.

*/
