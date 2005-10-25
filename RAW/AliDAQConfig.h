#ifndef ALIDAQCONFIG_H
#define ALIDAQCONFIG_H
/* Copyright(c) 1998-1999, ALICE Experiment at CERN, All rights reserved. *
 * See cxx source for full Copyright notice                               */

//////////////////////////////////////////////////////////////////////////////
//                                                                          //
// The header file contains all the information about the DAQ configuration.//
// This includes the number of DDLs and LDC per sub-detector readout system.//
// This header file is to be included inside DATE.                          //
//                                                                          //
//////////////////////////////////////////////////////////////////////////////

static const int kNDetectors = 18;

static const char* kDetectors[kNDetectors] = {"TPC", "ITSSPD", "ITSSDD", "ITSSSD", 
					      "TRD", "TOF", "PHOS", "RICH", 
					      "EMCAL", "MUON", "MUTR", "ZDC", 
					      "PMD", "START", "VZERO", "CRT",
					      "FMD","TRIG"};
static const int kDetectorDDLs[kNDetectors]   = {216, 20, 12, 16, 
						 18, 72, 20, 20, 
						 22, 20, 2, 1, 
						 6, 1, 1, 1,
						 3,1};
static const float kDetectorLDCs[kNDetectors] = {46, 2, 2, 1, 
						 4, 2, 1, 2, 
						 1, 2, 1, 1,
						 1, 0.5, 0.5, 1,
						 1,1};


#endif
