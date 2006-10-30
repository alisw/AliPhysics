#ifdef __CINT__
/* Copyright(c) 1998-1999, ALICE Experiment at CERN, All rights reserved. *
 * See cxx source for full Copyright notice                               */

/* $Id$ */

#pragma link off all globals;
#pragma link off all classes;
#pragma link off all functions;
//#pragma link C++ enum   Cluster_t;

//#pragma link C++ global gITSdisplay;  // global used by AliITSdisplay

// Standard ITS classes 
 
#pragma link C++ class  AliITSRawCluster+;
#pragma link C++ class  AliITSRawClusterSPD+;
#pragma link C++ class  AliITSRawClusterSDD+;
#pragma link C++ class  AliITSRawClusterSSD+;
#pragma link C++ class  AliITSClusterFinder+;
#pragma link C++ class  AliITSClusterFinderSPD+;
//#pragma link C++ class  AliITSClusterFinderSPDdubna+;
#pragma link C++ class  AliITSClusterFinderSDD+;
#pragma link C++ class  AliITSClusterFinderSSD+;
#pragma link C++ class  AliITSDetTypeRec+;

#pragma link C++ class  AliITSclusterSSD+;
#pragma link C++ class  AliITSpackageSSD+;
// Classes used for Tracking
//#pragma link C++ class  AliITSTrackV1+;
#pragma link C++ class  AliITSRad+;
#pragma link C++ class  AliITSIOTrack+;
//#pragma link C++ class  AliITSTrackerV1+;
#pragma link C++ class  AliITSRiemannFit+;

#pragma link C++ class AliITSclustererV2+;
#pragma link C++ class AliITStrackV2+;
#pragma link C++ class AliITStrackerV2+;
#pragma link C++ class AliITStrackMI+;
#pragma link C++ class AliITStrackerMI+;
//#pragma link C++ class AliITSRecV0Info+;

#pragma link C++ class  AliITSVertexer+;
#pragma link C++ class  AliITSVertexerIons+;
#pragma link C++ class  AliITSVertexerPPZ+;
#pragma link C++ class  AliITSVertexerZ+;
#pragma link C++ class  AliITSVertexer3D+;

// Classes for neural tracking
#pragma link C++ class AliITSNeuralPoint+;
#pragma link C++ class AliITSNeuralTrack+;
#pragma link C++ class AliITSNeuralTracker+;
#pragma link C++ class AliITStrackerANN+;
// Tasks
#pragma link C++ class AliITSreconstruction+;
#pragma link C++ class AliITSFindClustersV2+;
//#pragma link C++ class DisplayITSv11+;

#pragma link C++ class AliITSclusterTable+;
#pragma link C++ class AliITStrackerSA+;
#pragma link C++ class AliITStrackSA+;
#pragma link C++ class AliITSVertexerFast+;
#pragma link C++ class AliITSReconstructor+;
#pragma link C++ class AliITSClusterFinderV2+;
#pragma link C++ class AliITSClusterFinderV2SDD+;
#pragma link C++ class AliITSClusterFinderV2SPD+;
#pragma link C++ class AliITSClusterFinderV2SSD+;

// Classes for PID
#pragma link C++ class  AliITSPid+;
#pragma link C++ class  AliITStrackV2Pid+;
#pragma link C++ class AliITSPidParItem+;
#pragma link C++ class AliITSPident+;
#pragma link C++ class AliITSSteerPid+;
#pragma link C++ class AliITSpidESD+;
#pragma link C++ class AliITSpidESD1+;
#pragma link C++ class AliITSpidESD2+;
//beam test classes
#pragma link C++ class AliITSBeamTestDig+;
#pragma link C++ class AliITSBeamTestDigSPD+;
#pragma link C++ class AliITSBeamTestDigSDD+;
#pragma link C++ class AliITSBeamTestDigSSD+;
#pragma link C++ class AliITSBeamTestDigitizer+;
//multiplicity
#pragma link C++ class AliITSMultReconstructor+;
// SPD and SDD preprocessing
#pragma link C++ class AliITSBadChannelsAuxSPD+;
#pragma link C++ class AliITSBadChannelsSPD+;
#pragma link C++ class AliITSChannelSPD+;
#pragma link C++ class AliITSPreprocessorSPD+;
#pragma link C++ class AliITSPreprocessorSDD+;

#endif
