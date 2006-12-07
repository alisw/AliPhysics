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

/*
$Log$
Revision 1.4  2006/11/29 16:08:01  hristov
START becomes T0

Revision 1.3  2006/11/29 09:53:27  hristov
RICH becomes HMPID

Revision 1.2  2006/11/06 14:24:21  jgrosseo
reading of run parameters from the logbook
online offline naming conversion

Revision 1.1  2006/06/02 14:14:36  hristov
Separate library for CDB (Jan)

Revision 1.2  2006/03/07 07:52:34  hristov
New version (B.Yordanov)

Revision 1.3  2005/11/17 17:47:34  byordano
TList changed to TObjArray

Revision 1.2  2005/11/17 14:43:22  byordano
import to local CVS

Revision 1.1.1.1  2005/10/28 07:33:58  hristov
Initial import as subdirectory in AliRoot

Revision 1.1.1.1  2005/09/12 22:11:40  byordano
SHUTTLE package

Revision 1.2  2005/08/29 21:15:47  byordano
some docs added

*/

//
// abstract interface class to AliShuttle
// This class is implemented by AliTestShuttle for testing and
// by AliShuttle for the full setup
//

#include "AliShuttleInterface.h"
#include "AliLog.h"
#include <TSystem.h>

ClassImp(AliShuttleInterface)

const char* AliShuttleInterface::fkSystemNames[3] = { "DAQ", "DCS", "HLT" };

// names of the detectors preprocessors
const char* AliShuttleInterface::fgkDetName[kNDetectors] = {"SPD", "SDD", "SSD", "TPC", "TRD", "TOF",
       "PHS", "CPV", "HMP", "EMC", "MCH", "MTR", "FMD", "ZDC", "PMD", "T00", "V00", "GRP"};

// names of the detectors in OCDB
const char* AliShuttleInterface::fgkOfflineDetName[kNDetectors] = {"ITS", "ITS", "ITS", "TPC", "TRD", "TOF",
       "PHOS", "PHOS", "HMPID", "EMCAL", "MUON", "MUON", "FMD", "ZDC", "PMD", "T0", "VZERO", "GRP"};

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
const Int_t AliShuttleInterface::GetDetPos(const char* detName){
// Return detector position in the detector code array

	for(UInt_t iDet=0; iDet < kNDetectors; iDet++){
		if(!strcmp(fgkDetName[iDet], detName)) return iDet;
	}
	return -1;
}
