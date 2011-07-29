#ifdef __CINT__
/* Copyright(c) 1998-1999, ALICE Experiment at CERN, All rights reserved. *
 * See cxx source for full Copyright notice                               */

/* $Id: ITSsimLinkDef.h 36856 2009-11-16 16:17:07Z masera $ */

#pragma link off all globals;
#pragma link off all classes;
#pragma link off all functions;
//#pragma link C++ enum   Cluster_t;

//#pragma link C++ global gITSdisplay;  // global used by AliITSdisplay

// ITS upgrade classes
#pragma link C++ class  AliITSupgrade+;
#pragma link C++ class  AliITSupgradeDigitizer+;
#pragma link C++ class  AliITSUpg+;
#pragma link C++ class  AliITSvUpgrade+;
#pragma link C++ class  AliITSv11GeometryUpgrade+;

#endif
