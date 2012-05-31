// $Id$

//**************************************************************************
//* This file is property of and copyright by the ALICE HLT Project        * 
//* ALICE Experiment at CERN, All rights reserved.                         *
//*                                                                        *
//* Primary Authors: Matthias Richter <Matthias.Richter@ift.uib.no>        *
//*                  for The ALICE HLT Project.                            *
//*                                                                        *
//* Permission to use, copy, modify and distribute this software and its   *
//* documentation strictly for non-commercial purposes is hereby granted   *
//* without fee, provided that the above copyright notice appears in all   *
//* copies and that both the copyright notice and this permission notice   *
//* appear in the supporting documentation. The authors make no claims     *
//* about the suitability of this software for any purpose. It is          *
//* provided "as is" without express or implied warranty.                  *
//**************************************************************************

/// @file   AliHLTDAQ.cxx
/// @author Matthias Richter
/// @date   24.10.2008
/// @brief  Virtual Interface to the AliDAQ class.
///

#include "AliHLTDAQ.h"
#include "AliHLTLogging.h"
#include "AliHLTDataTypes.h"
#include "AliHLTMisc.h"
#include "TClass.h"
#include "TSystem.h"

/** ROOT macro for the implementation of ROOT specific class methods */
ClassImp(AliHLTDAQ)

AliHLTDAQ::AliHLTDAQ()
{
  // see header file for class documentation
  // or
  // refer to README to build package
  // or
  // visit http://web.ift.uib.no/~kjeks/doc/alice-hlt
}

AliHLTDAQ* AliHLTDAQ::fgpInstance=NULL;
const char* AliHLTDAQ::fgkImplName="AliHLTDAQInterfaceImplementation";
const char* AliHLTDAQ::fgkImplLibrary="libHLTrec.so";
const char* AliHLTDAQ::fgkOriginMapping[] = {
  kAliHLTDataOriginITSSPD,  
  kAliHLTDataOriginITSSDD,
  kAliHLTDataOriginITSSSD,
  kAliHLTDataOriginTPC,
  kAliHLTDataOriginTRD,
  kAliHLTDataOriginTOF,
  kAliHLTDataOriginHMPID,
  kAliHLTDataOriginPHOS,
  kAliHLTDataOriginCPV,
  kAliHLTDataOriginPMD,
  kAliHLTDataOriginMUON,
  "", // MUONTRG
  kAliHLTDataOriginFMD,
  kAliHLTDataOriginT0,
  kAliHLTDataOriginVZERO,
  kAliHLTDataOriginZDC,
  kAliHLTDataOriginACORDE,
  kAliHLTDataOriginTRG,
  kAliHLTDataOriginEMCAL,
  NULL
};

AliHLTDAQ::~AliHLTDAQ()
{
  // see header file for class documentation
}

Int_t       AliHLTDAQ::NumberOfDetectors()
{
  // see header file for class documentation
  if (!fgpInstance) GetInstance();
  if (fgpInstance) return fgpInstance->VirtNumberOfDetectors();
  return -1;
}

Int_t       AliHLTDAQ::DetectorID(const char *detectorName)
{
  // get the detector ID from the detector name
  // forwards to AliDAQ::DetectorName
  if (!fgpInstance) GetInstance();
  if (fgpInstance) return fgpInstance->VirtDetectorID(detectorName);
  return -1;
}

const char *AliHLTDAQ::DetectorName(Int_t detectorID)
{
  // get the detector name from the detector ID
  // forwards to AliDAQ::DetectorName
  if (!fgpInstance) GetInstance();
  if (fgpInstance) return fgpInstance->VirtDetectorName(detectorID);
  return NULL;
}

Int_t       AliHLTDAQ::DetectorIDFromHLTOrigin(const char dataorigin[kAliHLTComponentDataTypefOriginSize])
{
  // get the detector id from HLT origin
  for (int i=0; fgkOriginMapping[i]!=NULL; i++) {
    if (strncmp(fgkOriginMapping[i], dataorigin, kAliHLTComponentDataTypefOriginSize)==0) {
      return i;
    }
  }

  return -1;
}

const char *AliHLTDAQ::DetectorName(const char dataorigin[kAliHLTComponentDataTypefOriginSize])
{
  // get the detector name from HLT origin
  Int_t detectorID=DetectorIDFromHLTOrigin(dataorigin);
  if (detectorID<0) return NULL;

  return DetectorName(detectorID);
}

Int_t       AliHLTDAQ::DdlIDOffset(const char *detectorName)
{
  // get ID of the first DDL of a detector
  // forwards to AliDAQ::DdlIDOffset
  if (!fgpInstance) GetInstance();
  if (fgpInstance) return fgpInstance->VirtDdlIDOffset(detectorName);
  return -1;
}

Int_t       AliHLTDAQ::DdlIDOffset(Int_t detectorID)
{
  // get ID of the first DDL of a detector
  // forwards to AliDAQ::DdlIDOffset
  if (!fgpInstance) GetInstance();
  if (fgpInstance) return fgpInstance->VirtDdlIDOffset(detectorID);
  return -1;
}

const char *AliHLTDAQ::DetectorNameFromDdlID(Int_t ddlID, Int_t &ddlIndex)
{
  // get detector name from global DDL ID and index of DDL within the detector
  // forwards to AliDAQ::DetectorNameFromDdlID
  if (!fgpInstance) GetInstance();
  if (fgpInstance) return fgpInstance->VirtDetectorNameFromDdlID(ddlID, ddlIndex);
  return NULL;
}

Int_t       AliHLTDAQ::DetectorIDFromDdlID(Int_t ddlID, Int_t &ddlIndex)
{
  // get detector ID from global DDL ID and index of DDL within the detector
  // forwards to AliDAQ::DetectorNameFromDdlID
  if (!fgpInstance) GetInstance();
  if (fgpInstance) return fgpInstance->VirtDetectorIDFromDdlID(ddlID, ddlIndex);
  return -1;
}

Int_t       AliHLTDAQ::DdlID(const char *detectorName, Int_t ddlIndex)
{
  // get global DDL ID from detector name and index of DDL within the detector
  // forwards to AliDAQ::DdlID
  if (!fgpInstance) GetInstance();
  if (fgpInstance) return fgpInstance->VirtDdlID(detectorName, ddlIndex);
  return -1;
}

Int_t       AliHLTDAQ::DdlID(Int_t detectorID, Int_t ddlIndex)
{
  // get global DDL ID from detector ID and index of DDL within the detector
  // forwards to AliDAQ::DdlID
  if (!fgpInstance) GetInstance();
  if (fgpInstance) return fgpInstance->VirtDdlID(detectorID, ddlIndex);
  return -1;
}

const char *AliHLTDAQ::DdlFileName(const char *detectorName, Int_t ddlIndex)
{
  // see header file for class documentation
  if (!fgpInstance) GetInstance();
  if (fgpInstance) return fgpInstance->VirtDdlFileName(detectorName, ddlIndex);
  return NULL;
}

const char *AliHLTDAQ::DdlFileName(Int_t detectorID, Int_t ddlIndex)
{
  // see header file for class documentation
  if (!fgpInstance) GetInstance();
  if (fgpInstance) return fgpInstance->VirtDdlFileName(detectorID, ddlIndex);
  return NULL;
}

Int_t       AliHLTDAQ::NumberOfDdls(const char *detectorName)
{
  // see header file for class documentation
  if (!fgpInstance) GetInstance();
  if (fgpInstance) return fgpInstance->VirtNumberOfDdls(detectorName);
  return -1;
}

Int_t       AliHLTDAQ::NumberOfDdls(Int_t detectorID)
{
  // see header file for class documentation
  if (!fgpInstance) GetInstance();
  if (fgpInstance) return fgpInstance->VirtNumberOfDdls(detectorID);
  return -1;
}

const char *AliHLTDAQ::ListOfTriggeredDetectors(UInt_t detectorPattern)
{
  // see header file for class documentation
  if (!fgpInstance) GetInstance();
  if (fgpInstance) return fgpInstance->VirtListOfTriggeredDetectors(detectorPattern);
  return NULL;
}

UInt_t      AliHLTDAQ::DetectorPattern(const char *detectorList)
{
  // see header file for class documentation
  if (!fgpInstance) GetInstance();
  if (fgpInstance) return fgpInstance->VirtDetectorPattern(detectorList);
  return 0;
}

const char *AliHLTDAQ::OfflineModuleName(const char *detectorName)
{
  // see header file for class documentation
  if (!fgpInstance) GetInstance();
  if (fgpInstance) return fgpInstance->VirtOfflineModuleName(detectorName);
  return NULL;
}
const char *AliHLTDAQ::OfflineModuleName(Int_t detectorID)
{
  // see header file for class documentation
  if (!fgpInstance) GetInstance();
  if (fgpInstance) return fgpInstance->VirtOfflineModuleName(detectorID);
  return NULL;
}

const char *AliHLTDAQ::OnlineName(const char *detectorName)
{
  // see header file for class documentation
  if (!fgpInstance) GetInstance();
  if (fgpInstance) return fgpInstance->VirtOnlineName(detectorName);
  return NULL;
}
const char *AliHLTDAQ::OnlineName(Int_t detectorID)
{
  // see header file for class documentation
  if (!fgpInstance) GetInstance();
  if (fgpInstance) return fgpInstance->VirtOnlineName(detectorID);
  return NULL;
}

string AliHLTDAQ::HLTOrigin(const char *detectorName)
{
  // get HLT origin from detector name
  return HLTOrigin(DetectorID(detectorName));
}

string AliHLTDAQ::HLTOrigin(Int_t detectorID)
{
  // get HLT origin from detector ID
  string origin;
  if (detectorID>=0 && (unsigned)detectorID<sizeof(fgkOriginMapping)/sizeof(const char*)) {
    origin.append(fgkOriginMapping[detectorID], kAliHLTComponentDataTypefOriginSize);
  }
  return origin;
}

string AliHLTDAQ::HLTSpecificationFromDdlID(Int_t ddlID)
{
  // get the HLT data specification from equipment id
  string result = "";
  Int_t ddlIndex=0;
  Int_t detectorID = DetectorIDFromDdlID(ddlID, ddlIndex);
  if (detectorID < 0)
    return result;
  Int_t TPCID = DetectorID("TPC");
  const int strtmplength=11;
  char strtmp[strtmplength];
  memset(strtmp, 0, strtmplength);
  if (detectorID == TPCID) {
    int partition;
    int slice;
    if (ddlID < 840) {
      partition = ddlID % 2;
      slice = (ddlID - 768) / 2;
    } else {
      partition = (ddlID % 4) + 2;
      slice = (ddlID - 840) / 4;
    }
    snprintf(strtmp, strtmplength, "0x%02x%02x%02x%02x", slice, slice, partition, partition);
    result = strtmp;
  }
  else if (detectorID == DetectorID("TOF")) {
    AliHLTLogging log;
    log.Logging(kHLTLogWarning, "AliHLTDAQ::HLTSpecificationFromDdlID", "HLT Analysis", "Mapping of data specification not implemented for TOF");
  }
  else { // default
    snprintf(strtmp, strtmplength, "0x%08x", 0x1 << ddlIndex);
    result = strtmp;
  }
  return result;
}

Int_t AliHLTDAQ::DdlIDFromHLTBlockData(const char dataorigin[kAliHLTComponentDataTypefOriginSize], UInt_t specification)
{
  // get the DDL ID (global equipment ID) from HLT origin and data specification
  Int_t detectorID=DetectorIDFromHLTOrigin(dataorigin);
  if (detectorID<0) return -1;
  Int_t ddlID=DdlIDOffset(detectorID);
  if (ddlID<0) return -1;

  if (detectorID == DetectorID("TPC")) {
    int minPartition= specification     &0xff;
    int maxPartition=(specification>>8) &0xff;
    int minSlice    =(specification>>16)&0xff;
    int maxSlice    =(specification>>24)&0xff;
    if (minPartition<0 || minPartition>5 ||
	maxPartition<0 || maxPartition>5 ||
	minSlice<0 || minSlice>35 ||
	maxSlice<0 || maxSlice>35) {
      AliHLTLogging log;
      log.Logging(kHLTLogError, "AliHLTDAQ::DdlID", "HLT Analysis", "invalid data specification 0x%08x", specification);
      return -1;
    }
    else if (minPartition!=maxPartition || 
	     minSlice!=maxSlice) {
      AliHLTLogging log;
      log.Logging(kHLTLogError, "AliHLTDAQ::DdlID", "HLT Analysis", "ambiguous data specification 0x%08x", specification);
      return -1;
    }
    if (minPartition<2) ddlID+=2*minSlice+minPartition;        // inner sector
    else                ddlID+=4*minSlice+(minPartition-2)+72; // outer sector
  }
  else if (detectorID == DetectorID("TOF")) {
    AliHLTLogging log;
    log.Logging(kHLTLogWarning, "AliHLTDAQ::DdlID", "HLT Analysis", "Mapping of data specification not implemented for TOF");
    return -1;
  }
  else { // default
    for (int i=0; i<32; i++) {
      if ((specification&(0x1<<i))==0) continue; // if bit not set
      if ((specification&(0xffffffff^(0x1<<i)))!=0) { // check if other bits set
	AliHLTLogging log;
	log.Logging(kHLTLogError, "AliHLTDAQ::DdlID", "HLT Analysis", "ambiguous data specification 0x%08x", specification);
	return -1;
      }
      ddlID+=i;
      break;
    }
  }

  return ddlID;
}

AliHLTDAQ* AliHLTDAQ::GetInstance()
{
  // see header file for class documentation
  if (!fgpInstance) {
    fgpInstance=AliHLTMisc::LoadInstance((AliHLTDAQ*)NULL, fgkImplName, fgkImplLibrary);
    if (!fgpInstance) {
      AliHLTLogging log;
      log.Logging(kHLTLogError, "AliHLTDAQ::Instance", "HLT Analysis", "failed to load AliHLTDAQ implementation");
    }
  }
  return fgpInstance;
}
