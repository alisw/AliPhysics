/* Copyright(c) 1998-1999, ALICE Experiment at CERN, All rights reserved. *
 * See cxx source for full Copyright notice                               */

// $Id$

/// \file MUONrecLinkDef.h
/// \brief The CINT link definitions for \ref rec 

#ifdef __CINT__
#pragma link off all globals;
#pragma link off all classes;
#pragma link off all functions;

// reconstruction
#pragma link C++ class AliMUONReconstructor+;
#pragma link C++ class AliMUONVTrackReconstructor+; 
#pragma link C++ class AliMUONTrackReconstructor+; 
#pragma link C++ class AliMUONTrackReconstructorK+; 
#pragma link C++ class AliMUONTracker+;
#pragma link C++ class AliMUONTrack+; 
#pragma link C++ class AliMUONTrackParam+; 
#pragma link C++ class AliMUONTrackExtrap+; 
#pragma link C++ class AliMUONTriggerTrack+; 
#pragma link C++ class AliMUONRecoTrack+; 
#pragma link C++ class AliMUONAlignment+;
#pragma link C++ class AliMUONVClusterFinder+;
#pragma link C++ class AliMUONPad+;
#pragma link C++ class AliMUONCluster+;
#pragma link C++ class AliMUONPreClusterFinder+;
#pragma link C++ class AliMUONPreClusterFinderV2+;
#pragma link C++ class AliMUONPreClusterFinderV3+;
#pragma link C++ class AliMUONClusterFinderSimpleFit+;
#pragma link C++ class AliMUONClusterFinderCOG+;
#pragma link C++ class AliMUONClusterFinderPeakCOG+;
#pragma link C++ class AliMUONClusterFinderPeakFit+;
#pragma link C++ class AliMUONClusterFinderMLEM+;
#pragma link C++ class AliMUONClusterSplitterMLEM+;
#pragma link C++ class AliMUONTrackHitPattern+;
#pragma link C++ class AliMUONRefitter+;

#pragma link C++ class AliMUONVClusterStore+;
#pragma link C++ class AliMUONClusterStoreV1+;
#pragma link C++ class AliMUONClusterStoreV2+;
#pragma link C++ class AliMUONClusterStoreV2Iterator+;
#pragma link C++ class AliMUONVTrackStore+;
#pragma link C++ class AliMUONTrackStoreV1+;
#pragma link C++ class AliMUONVTriggerTrackStore+;
#pragma link C++ class AliMUONTriggerTrackStoreV1+;

#pragma link C++ class AliMUONRecoParam+;

#pragma link C++ class AliMUONVClusterServer+;
#pragma link C++ class AliMUONSimpleClusterServer+;
#pragma link C++ class AliMUONLegacyClusterServer+;
#pragma link C++ class AliMUONTriggerTrackToTrackerClusters+;

#pragma link C++ class AliMUONESDInterface+;

// calibration
#pragma link C++ class AliMUONDigitCalibrator+;
#pragma link C++ class AliMUONPadStatusMaker+;
#pragma link C++ class AliMUONPadStatusMapMaker+;

// QA
#pragma link C++ class AliMUONQADataMakerRec+;
#pragma link C++ class AliMUONQAChecker+;

#pragma link C++ class AliMUONVTrackerDataMaker+;
#pragma link C++ class AliMUONTrackerDataMaker+;

#endif
