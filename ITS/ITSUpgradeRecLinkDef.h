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
 
#pragma link C++ class AliITSlayerUpgrade+;
#pragma link C++ class AliITStrackerUpgrade+;
//#pragma link C++ class AliITSReconstructor+;
#pragma link C++ class AliITSUpgradeReconstructor+;
#pragma link C++ class AliITSUpgradeClusterList+;
#pragma link C++ class AliITSUpgradeClusterListNode+;
#pragma link C++ class AliITSUpgradeClusterFinder+;

#endif
