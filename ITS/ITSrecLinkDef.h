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
 
#pragma link C++ class  AliITSClusterFinder+;
#pragma link C++ class  AliITSClusterParam+;
#pragma link C++ class  AliITSDetTypeRec+;

#pragma link C++ class  AliITSclusterSSD+;
#pragma link C++ class  AliITSpackageSSD+;
// Classes used for Tracking
//#pragma link C++ class  AliITSTrackV1+;
#pragma link C++ class  AliITSRad+;
#pragma link C++ class  AliITSIOTrack+;
//#pragma link C++ class  AliITSTrackerV1+;

#pragma link C++ class AliITSclustererV2+;
#pragma link C++ class AliITStrackV2+;
#pragma link C++ class AliITStrackerV2+;
#pragma link C++ class AliITStrackMI+;
#pragma link C++ class AliITStrackerMI+;
//#pragma link C++ class AliITSRecV0Info+;

#pragma link C++ class  AliITSVertexer+;
#pragma link C++ class  AliITSVertexerIons+;
#pragma link C++ class  AliITSVertexerCosmics+;
#pragma link C++ class  AliITSVertexerZ+;
#pragma link C++ class  AliITSVertexer3D+;
#pragma link C++ class  AliITSVertexer3DTapan+;
#pragma link C++ class  AliITSTracklPairs+;
#pragma link C++ class  AliITSSortTrkl+;
#pragma link C++ class AliITSVertexerFast+;
#pragma link C++ class AliITSVertexerFixed+;
#pragma link C++ class  AliITSMeanVertexer+;
#pragma link C++ class  AliITSZPoint+;

// Tasks
#pragma link C++ class AliITSreconstruction+;
//#pragma link C++ class DisplayITSv11+;

#pragma link C++ class AliITSclusterTable+;
#pragma link C++ class AliITStrackerSA+;
#pragma link C++ class AliITStrackSA+;
#pragma link C++ class AliITSReconstructor+;
#pragma link C++ class AliITSRecoParam+;
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
//multiplicity with tracklets
#pragma link C++ class AliITSTrackleterSPDEff+;
#pragma link C++ class AliITSMultReconstructor+;

// SPD, SDD and SSD preprocessing
#pragma link C++ class AliITSBadChannelsAuxSPD+;
#pragma link C++ class AliITSBadChannelsSPD+;
#pragma link C++ class AliITSChannelSPD+;
#pragma link C++ class AliITSPreprocessorSPD+;
#pragma link C++ class AliITSOnlineCalibrationSPD+;
#pragma link C++ class AliITSOnlineCalibrationSPDhandler+;
#pragma link C++ class AliITSOnlineSPDHitArray+;
#pragma link C++ class AliITSOnlineSPDHitEvent+;
#pragma link C++ class AliITSOnlineSPDscanAnalyzer+;
#pragma link C++ class AliITSOnlineSPDscan+;
#pragma link C++ class AliITSOnlineSPDscanInfo+;
#pragma link C++ class AliITSOnlineSPDscanInfoMeanTh+;
#pragma link C++ class AliITSOnlineSPDscanInfoMultiple+;
#pragma link C++ class AliITSOnlineSPDscanMeanTh+;
#pragma link C++ class AliITSOnlineSPDscanMultiple+;
#pragma link C++ class AliITSOnlineSPDscanSingle+;
#pragma link C++ class AliITSOnlineSPDphys+;
#pragma link C++ class AliITSOnlineSPDphysAnalyzer+;
#pragma link C++ class AliITSOnlineSPDphysInfo+;
#pragma link C++ class AliITSOnlineSPDfoChipConfig+;
#pragma link C++ class AliITSOnlineSPDfoChip+;
#pragma link C++ class AliITSOnlineSPDfo+;
#pragma link C++ class AliITSOnlineSPDfoInfo+;
#pragma link C++ class AliITSOnlineSPDfoAnalyzer+;
#pragma link C++ class AliITSPreprocessorSDD+;
#pragma link C++ class AliITSOnlineSDD+;
#pragma link C++ class AliITSOnlineSDDBase+;
#pragma link C++ class AliITSOnlineSDDTP+;
#pragma link C++ class AliITSOnlineSDDCMN+;
#pragma link C++ class AliITSOnlineSDDInjectors+;
#pragma link C++ class AliITSPreprocessorSSD+;

// Classes for alignment
#pragma link C++ class AliITSAlignMille+;
#pragma link C++ class AliITSAlignMille2+;
#pragma link C++ class AliITSAlignMilleModule+;
#pragma link C++ class AliITSAlignMilleData+;
#pragma link C++ class AliITSAlignMille2Module+;
#pragma link C++ class AliITSResidualsAnalysis+;
#pragma link C++ class AliITSRealignTracks+;
// Classes for QA
#pragma link C++ class AliITSQAChecker+;
#pragma link C++ class AliITSQADataMakerRec+;
#pragma link C++ class AliITSQASPDDataMakerRec+;
#pragma link C++ class AliITSQASDDDataMakerRec+;
#pragma link C++ class AliITSQASSDDataMakerRec+;
#pragma link C++ class AliITSQASPDChecker+;
#pragma link C++ class AliITSQASDDChecker+;
#pragma link C++ class AliITSQASSDChecker+;

#endif
