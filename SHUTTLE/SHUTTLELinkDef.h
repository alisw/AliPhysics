#ifdef __CINT__
/* Copyright(c) 1998-1999, ALICE Experiment at CERN, All rights reserved. *
 * See cxx source for full Copyright notice                               */

/* $Id$ */

#pragma link off all globals;
#pragma link off all classes;
#pragma link off all functions;

#pragma link C++ global gAlice;
#pragma link C++ global gMC;
 
// AliDCSClient classes ...
#pragma link C++ class AliDCSClient+;
#pragma link C++ class AliDCSMessage+;

// Shuttle classes ...
#pragma link C++ class  AliShuttleConfig;
#pragma link C++ class  AliShuttleConfig::AliShuttleConfigHolder;
#pragma link C++ class  AliShuttle;
#pragma link C++ class  AliShuttleTrigger;
#pragma link C++ class  TerminateSignalHandler;
#pragma link C++ class  AliShuttleStatus;
#pragma link C++ class  AliShuttleLogbookEntry;

#endif
