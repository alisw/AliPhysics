#ifdef __CINT__
/* Copyright(c) 1998-1999, ALICE Experiment at CERN, All rights reserved. *
 * See cxx source for full Copyright notice                               */

/* $Id: ITSrecLinkDef.h 38329 2010-01-17 19:17:24Z hristov $ */

#pragma link off all globals;
#pragma link off all classes;
#pragma link off all functions;
//#pragma link C++ enum   Cluster_t;

//#pragma link C++ global gITSdisplay;  // global used by AliITSdisplay

// ITS upgrade classes 
 
#pragma link C++ class AliITSURecoParam+;
#pragma link C++ class AliITSUReconstructor+;
#pragma link C++ class AliITSUClusterizer+;
#pragma link C++ class AliITSUClusterLines+;
#pragma link C++ class AliITSURecoSens+;
#pragma link C++ class AliITSURecoLayer+;
#pragma link C++ class AliITSURecoDet+;
#pragma link C++ class AliITSUClusterPix+;
#pragma link C++ class AliITSUSeed+;
#pragma link C++ class AliITSUTrackerGlo+;
#pragma link C++ class AliITSUCATracker+;
#pragma link C++ class AliITSUTrackerCooked+;
#pragma link C++ class AliITSUTrackCooked+;
#pragma link C++ class AliITSUTrackCond+;
#pragma link C++ class AliITSUTrackHyp+;
#pragma link C++ class AliITSUMatLUT+;
#pragma link C++ class AliITSUVertexer+;

//

// old v0
/*#pragma link C++ class AliITSlayerUpgrade+;
#pragma link C++ class AliITStrackerUpgrade+;
#pragma link C++ class AliITStrackU+;
#pragma link C++ class AliITStrackerU+;
#pragma link C++ class AliITSUpgradeReconstructor+;
#pragma link C++ class AliITSUpgradeClusterList+;
#pragma link C++ class AliITSUpgradeClusterListNode+;
#pragma link C++ class AliITSUPixelChip+;
#pragma link C++ class AliITSUpgradeClusterFinder+;
*/

#endif
