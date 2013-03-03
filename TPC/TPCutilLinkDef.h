#ifdef __CINT__
/* Copyright(c) 1998-1999, ALICE Experiment at CERN, All rights reserved. *
 * See cxx source for full Copyright notice                               */

/* $Id: TPCutilLinkDef.h  $ */

#pragma link off all globals;
#pragma link off all classes;
#pragma link off all functions;

#pragma link C++ class AliTPCAltroEmulator+;   //UTIL Altro emulator -- Used 2010 for checks for MAF and TCF in Altro does proper job
#pragma link C++ class AliTPCAlign;            //UTIL Transform survey measurements to TPC alignment
#pragma link C++ class AliTPCCalibViewerGUIAlarms+;            //UTIL Used in Expert AMORE : Visualize Outlier histograms based on reference data - useful for use after LS1
#pragma link C++ class AliTPCCalibQAChecker+;                  //UTIL Used in Expert AMORE : Does the checking for above - useful for use after LS1
#pragma link C++ class AliTPCkalmanTime+;      //UTIL  Estimate coefficient for Time dependent calibration (for studies)
#pragma link C++ class AliTPCGenDBTemp+;       //UTIL  Used to generate OCDB entries by hand (Haavard) --- add documentation
#pragma link C++ class AliTPCGenDBConf+;       //UTIL  Used to generate OCDB entries by hand (Haavard) --- add documentation

#endif
