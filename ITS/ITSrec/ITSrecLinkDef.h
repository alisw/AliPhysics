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

// Classes used for Tracking
//#pragma link C++ class  AliITSTrackV1+;
#pragma link C++ class  AliITSRad+;
#pragma link C++ class  AliITSIOTrack+;
//#pragma link C++ class  AliITSTrackerV1+;

#pragma link C++ class AliITSclustererV2+;
//#pragma link C++ class AliITSRecV0Info+;

#pragma link C++ class  AliITSVertexer3DTapan+;
#pragma link C++ class  AliITSMeanVertexer+;

// Tasks
#pragma link C++ class AliITSreconstruction+;
//#pragma link C++ class DisplayITSv11+;

#pragma link C++ class AliITSCorrectSDDPoints+;

// Classes for PID
#pragma link C++ class AliITSdEdxAnalyzer+;
//multiplicity with tracklets

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
#pragma link C++ class AliITSAlignMille2Module+;
#pragma link C++ class AliITSAlignMille2Constraint+;
#pragma link C++ class AliITSAlignMille2ConstrArray+;
#pragma link C++ class AliITSAlignMilleData+;
#pragma link C++ class AliITSTPArrayFit+;
#pragma link C++ class AliITSRealignTracks+;
#pragma link C++ class AliITSSumTP+;
// Classes for QA
#pragma link C++ class AliITSQASSDRefData+;

#endif
