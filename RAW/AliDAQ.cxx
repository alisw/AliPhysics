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

//////////////////////////////////////////////////////////////////////////////
//                                                                          //
// The AliDAQ class is responsible for handling all the information about   //
// Data Acquisition configuration. It defines the detector indexing,        //
// the number of DDLs and LDCs per detector.                                //
// The number of LDCs per detector is used only in the simulation in order  //
// to define the configuration of the dateStream application. Therefore the //
// numbers in the corresponding array can be changed without affecting the  //
// rest of the aliroot code.                                                //
// The equipment ID (DDL ID) is an integer (32-bit) number defined as:      //
// Equipment ID = (detectorID << 8) + DDLIndex                              //
// where the detectorID is given by fgkDetectorName array and DDLIndex is   //
// the index of the corresponding DDL inside the detector partition.        //
// Due to DAQ/HLT limitations, the ddl indexes should be consequtive, or    //
// at least without big gaps in between.                                    //
// The sub-detector code use only this class in the simulation and reading  //
// of the raw data.                                                         //
//                                                                          //
// cvetan.cheshkov@cern.ch  2006/06/09                                      //
//                                                                          //
//////////////////////////////////////////////////////////////////////////////

#include <TClass.h>
#include <TString.h>

#include "AliDAQ.h"
#include "AliLog.h"

ClassImp(AliDAQ)

const char* AliDAQ::fgkDetectorName[AliDAQ::kNDetectors] = {
  "ITSSPD",
  "ITSSDD",
  "ITSSSD",
  "TPC",
  "TRD",
  "TOF",
  "HMPID",  // Name to be changed to HMPID
  "PHOS",
  "CPV",
  "PMD",
  "MUONTRK",
  "MUONTRG",
  "FMD",
  "START", // Name to be changed to T0
  "VZERO", // Name to be changed to V0 ?
  "ZDC",
  "CRT",   // Name to be changed to ACCORDE
  "TRG",
  "EMCAL",
  "HLT"
};

Int_t AliDAQ::fgkNumberOfDdls[AliDAQ::kNDetectors] = {
  20,
  24,
  16,
  216,
  18,
  72,
  20,
  20,
  10,
  6,
  20,
  2,
  3,
  1,
  1,
  1,
  1,
  1,
  24,
  10
};

Float_t AliDAQ::fgkNumberOfLdcs[AliDAQ::kNDetectors] = {
  4,
  4,
  4,
  36,
  3,
  12,
  4,
  4,
  2,
  1,
  4,
  1,
  1,
  0.5,
  0.5,
  1,
  1,
  1,
  4,
  5
};

AliDAQ::AliDAQ(const AliDAQ& source) :
  TObject(source)
{
  // Copy constructor
  // Nothing to be done
}

AliDAQ& AliDAQ::operator = (const AliDAQ& /* source */)
{
  // Assignment operator
  // Nothing to be done
  return *this;
}

Int_t AliDAQ::DetectorID(const char *detectorName)
{
  // Return the detector index
  // corresponding to a given
  // detector name
  TString detStr = detectorName;

  Int_t iDet;
  for(iDet = 0; iDet < kNDetectors; iDet++) {
    if (detStr.CompareTo(fgkDetectorName[iDet],TString::kIgnoreCase) == 0)
      break;
  }
  if (iDet == kNDetectors) {
    AliErrorClass(Form("Invalid detector name: %s !",detectorName));
    return -1;
  }
  return iDet;
}

const char *AliDAQ::DetectorName(Int_t detectorID)
{
  // Returns the name of particular
  // detector identified by its index
  if (detectorID < 0 || detectorID >= kNDetectors) {
    AliErrorClass(Form("Invalid detector index: %d (%d -> %d) !",detectorID,0,kNDetectors-1));
    return "";
  }
  return fgkDetectorName[detectorID];
}

Int_t AliDAQ::DdlIDOffset(const char *detectorName)
{
  // Returns the DDL ID offset
  // for a given detector
  Int_t detectorID = DetectorID(detectorName);
  if (detectorID < 0)
    return -1;
  
  return DdlIDOffset(detectorID);
}

Int_t AliDAQ::DdlIDOffset(Int_t detectorID)
{
  // Returns the DDL ID offset
  // for a given detector identified
  // by its index
  if (detectorID < 0 || detectorID >= kNDetectors) {
    AliErrorClass(Form("Invalid detector index: %d (%d -> %d) !",detectorID,0,kNDetectors-1));
    return -1;
  }
  return (detectorID << 8);
}

const char *AliDAQ::DetectorNameFromDdlID(Int_t ddlID,Int_t &ddlIndex)
{
  // Returns the detector name for
  // a given DDL ID
  ddlIndex = -1;
  Int_t detectorID = DetectorIDFromDdlID(ddlID,ddlIndex);
  if (detectorID < 0)
    return "";

  return DetectorName(detectorID);
}

Int_t AliDAQ::DetectorIDFromDdlID(Int_t ddlID,Int_t &ddlIndex)
{
  // Returns the detector ID and
  // the ddl index within the
  // detector range for
  // a given input DDL ID
  Int_t detectorID = ddlID >> 8;
  if (detectorID < 0 || detectorID >= kNDetectors) {
    AliErrorClass(Form("Invalid detector index: %d (%d -> %d) !",detectorID,0,kNDetectors-1));
    return -1;
  }
  ddlIndex = ddlID & 0xFF;
  if (ddlIndex >= fgkNumberOfDdls[detectorID]) {
    AliErrorClass(Form("Invalid DDL index %d (%d -> %d) for detector %d",
		       ddlIndex,0,fgkNumberOfDdls[detectorID],detectorID));
    ddlIndex = -1;
    return -1;
  }
  return detectorID;
}

Int_t AliDAQ::DdlID(const char *detectorName, Int_t ddlIndex)
{
  // Returns the DDL ID starting from
  // the detector name and the DDL
  // index inside the detector
  Int_t detectorID = DetectorID(detectorName);
  if (detectorID < 0)
    return -1;

  return DdlID(detectorID,ddlIndex);
}

Int_t AliDAQ::DdlID(Int_t detectorID, Int_t ddlIndex)
{
  // Returns the DDL ID starting from
  // the detector ID and the DDL
  // index inside the detector
  Int_t ddlID = DdlIDOffset(detectorID);
  if (ddlID < 0)
    return -1;
  
  if (ddlIndex >= fgkNumberOfDdls[detectorID]) {
    AliErrorClass(Form("Invalid DDL index %d (%d -> %d) for detector %d",
		       ddlIndex,0,fgkNumberOfDdls[detectorID],detectorID));
    return -1;
  }

  ddlID += ddlIndex;
  return ddlID;
}

const char *AliDAQ::DdlFileName(const char *detectorName, Int_t ddlIndex)
{
  // Returns the DDL file name
  // (used in the simulation) starting from
  // the detector name and the DDL
  // index inside the detector
  Int_t detectorID = DetectorID(detectorName);
  if (detectorID < 0)
    return "";

  return DdlFileName(detectorID,ddlIndex);
}

const char *AliDAQ::DdlFileName(Int_t detectorID, Int_t ddlIndex)
{
  // Returns the DDL file name
  // (used in the simulation) starting from
  // the detector ID and the DDL
  // index inside the detector
  Int_t ddlID = DdlIDOffset(detectorID);
  if (ddlID < 0)
    return "";
  
  if (ddlIndex >= fgkNumberOfDdls[detectorID]) {
    AliErrorClass(Form("Invalid DDL index %d (%d -> %d) for detector %d",
		       ddlIndex,0,fgkNumberOfDdls[detectorID],detectorID));
    return "";
  }

  ddlID += ddlIndex;
  TString fileName = DetectorName(detectorID);
  fileName += "_";
  fileName += ddlID;
  fileName += ".ddl";
  return fileName.Data();
}

Int_t AliDAQ::NumberOfDdls(const char *detectorName)
{
  // Returns the number of DDLs for
  // a given detector
  Int_t detectorID = DetectorID(detectorName);
  if (detectorID < 0)
    return -1;

  return NumberOfDdls(detectorID);
}

Int_t AliDAQ::NumberOfDdls(Int_t detectorID)
{
  // Returns the number of DDLs for
  // a given detector
  if (detectorID < 0 || detectorID >= kNDetectors) {
    AliErrorClass(Form("Invalid detector index: %d (%d -> %d) !",detectorID,0,kNDetectors-1));
    return -1;
  }

  return fgkNumberOfDdls[detectorID];
}

Float_t AliDAQ::NumberOfLdcs(const char *detectorName)
{
  // Returns the number of DDLs for
  // a given detector
  Int_t detectorID = DetectorID(detectorName);
  if (detectorID < 0)
    return -1;

  return NumberOfLdcs(detectorID);
}

Float_t AliDAQ::NumberOfLdcs(Int_t detectorID)
{
  // Returns the number of DDLs for
  // a given detector
  if (detectorID < 0 || detectorID >= kNDetectors) {
    AliErrorClass(Form("Invalid detector index: %d (%d -> %d) !",detectorID,0,kNDetectors-1));
    return -1;
  }

  return fgkNumberOfLdcs[detectorID];
}

void AliDAQ::PrintConfig()
{
  // Print the DAQ configuration
  // for all the detectors
  printf("====================================================================\n"
	 "|                ALICE Data Acquisition Configuration              |\n"
	 "====================================================================\n"
	 "| Detector ID | Detector Name | DDL Offset | # of DDLs | # of LDCs |\n"
	 "====================================================================\n");
  for(Int_t iDet = 0; iDet < kNDetectors; iDet++) {
    printf("|%11d  |%13s  |%10d  |%9d  |%9.1f  |\n",
	   iDet,DetectorName(iDet),DdlIDOffset(iDet),NumberOfDdls(iDet),NumberOfLdcs(iDet));
  }
  printf("====================================================================\n");

}
