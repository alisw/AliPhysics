/**************************************************************************
 * Copyright(c) 1998-1999, ALICE Experiment at CERN, All rights reserved. *
 *                                                                        *
 * Author: The ALICE Off-line Project.                                    *
 * Contributors are mentioned in the code where appropriate.              *
 *                                                                        *
 * Permission to use, copy, modify and distribute this software and its   *
 * documentation strictly for non-commercial purposes is hereby granted   *
 * without fee, provided that the above copyright notice appears in all   *
 * copies and that both the copyright notice and this permission notice   *
 * appear in the supporting documentation. The authors make no claims     *
 * about the suitability of this software for any purpose. It is          *
 * provided "as is" without express or implied warranty.                  *
 **************************************************************************/

/* $Id$ */

//
// abstract interface class to AliShuttle
// This class is implemented by AliTestShuttle for testing and
// by AliShuttle for the full setup
//

#include "AliShuttleInterface.h"
#include "AliLog.h"
#include <TClass.h>
#include <TSystem.h>

ClassImp(AliShuttleInterface)

	const char* AliShuttleInterface::fkSystemNames[4] = { "DAQ", "DCS", "HLT", "DQM" };

// names of the detectors preprocessors
const char* AliShuttleInterface::fgkDetName[kNDetectors] = {"SPD", "SDD", "SSD", "TPC", "TRD", "TOF",
       "PHS", "CPV", "HMP", "EMC", "MCH", "MTR", "FMD", "ZDC", "PMD", "T00", "V00", "GRP", "HLT", "ACO", "TRI"
// #ifdef MFT_UPGRADE
// 							    , "MFT" 
// #endif 
							    , "MFT"     // AU
       };

// names of the detectors in OCDB
const char* AliShuttleInterface::fgkOfflineDetName[kNDetectors] = {"ITS", "ITS", "ITS", "TPC", "TRD", "TOF",
       "PHOS", "PHOS", "HMPID", "EMCAL", "MUON", "MUON", "FMD", "ZDC", "PMD", "T0", "VZERO", "GRP", "HLT", "ACORDE", "TRIGGER"
// #ifdef MFT_UPGRADE
// 	   , "MFT" 
// #endif 
								   , "MFT"     // AU
	   };

TString AliShuttleInterface::fgkMainCDB("alien://folder=ShuttleCDB");
TString AliShuttleInterface::fgkLocalCDB("local://LocalShuttleCDB");
TString AliShuttleInterface::fgkMainRefStorage("alien://folder=ShuttleReference");
TString AliShuttleInterface::fgkLocalRefStorage("local://LocalReferenceStorage");

TString AliShuttleInterface::fgkShuttleTempDir("/tmp");
TString AliShuttleInterface::fgkShuttleLogDir("/tmp/log");

//______________________________________________________________________________________________
const char* AliShuttleInterface::GetOfflineDetName(const char* detName){
// Return "offline" detector name

	Int_t detPos = GetDetPos(detName);
	if(detPos < 0) {
		AliErrorClass(Form("Unknown detector: %s",detName));
		return 0;
	}

	return fgkOfflineDetName[detPos];
}

//______________________________________________________________________________________________
const char* AliShuttleInterface::GetDetName(UInt_t detPos){
// Return detector code

	if(detPos >= kNDetectors) {
		AliErrorClass(Form("Parameter out of bound: %d", detPos));
		return 0;
	}

	return fgkDetName[detPos];
}

//______________________________________________________________________________________________
Int_t AliShuttleInterface::GetDetPos(const char* detName){
// Return detector position in the detector code array

	for(UInt_t iDet=0; iDet < kNDetectors; iDet++){
		if(!strcmp(fgkDetName[iDet], detName)) return iDet;
	}
	return -1;
}
