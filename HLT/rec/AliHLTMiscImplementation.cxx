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

/// @file   AliHLTMisc.cxx
/// @author Matthias Richter
/// @date   
/// @brief  Miscellaneous methods for the HLT AliRoot integration

#include "AliHLTMiscImplementation.h"
#include "AliHLTLogging.h"
#include "AliCDBManager.h"
#include "AliCDBStorage.h"
#include "AliCDBEntry.h"
#include "AliGRPManager.h"
#include "AliRawReader.h"
#include "TGeoGlobalMagField.h"

/** ROOT macro for the implementation of ROOT specific class methods */
ClassImp(AliHLTMiscImplementation);

AliHLTMiscImplementation::AliHLTMiscImplementation()
{
}

AliHLTMiscImplementation::~AliHLTMiscImplementation()
{
  // see header file for function documentation
}

int AliHLTMiscImplementation::InitCDB(const char* cdbpath)
{
  // see header file for function documentation
  int iResult=0;
  AliCDBManager* pCDB = AliCDBManager::Instance();
  AliHLTLogging log;
  if (!pCDB) {
    log.Logging(kHLTLogError, "InitCDB", "CDB handling", "Could not get CDB instance");
  } else {
    if (cdbpath && cdbpath[0]!=0) {
      TString uri=cdbpath;
      if (!uri.BeginsWith("local://") &&
	  !uri.Contains("://")) {
	// assume a local path if no uri specifier is found
	uri="local://";
	uri+=cdbpath;
      }
      pCDB->SetDefaultStorage(uri);
      log.Logging(kHLTLogDebug, "InitCDB", "CDB handling", "CDB instance 0x%x", pCDB);
    } else if (!pCDB->IsDefaultStorageSet()) {
      const char* cdbUri="local://$ALICE_ROOT/OCDB";
      pCDB->SetDefaultStorage(cdbUri);
      pCDB->SetRun(0);
      log.Logging(kHLTLogInfo, "InitCDB", "CDB handling", "set default URI: %s", cdbUri);
    }
  }
  return iResult;
}

int AliHLTMiscImplementation::SetCDBRunNo(int runNo)
{
  // see header file for function documentation
  int iResult=0;
  AliCDBManager* pCDB = AliCDBManager::Instance();
  AliHLTLogging log;
  if (!pCDB) {
    log.Logging(kHLTLogError, "SetCDBRunNo", "CDB handling", "Could not get CDB instance");
  } else {
    pCDB->SetRun(runNo);
  }
  return iResult;
}

AliCDBEntry* AliHLTMiscImplementation::LoadOCDBEntry(const char* path, int runNo, int version, int subVersion)
{
  // see header file for function documentation
  AliCDBStorage* store = AliCDBManager::Instance()->GetDefaultStorage();
  if (!store) {
    return NULL;
  }
  if (version<0) version = store->GetLatestVersion(path, runNo);
  if (subVersion<0) subVersion = store->GetLatestSubVersion(path, runNo, version);
  return AliCDBManager::Instance()->Get(path, runNo, version, subVersion);
}

TObject* AliHLTMiscImplementation::ExtractObject(AliCDBEntry* entry)
{
  // see header file for function documentation
  if (!entry) return NULL;
  return entry->GetObject();
}

int AliHLTMiscImplementation::InitMagneticField() const
{
  // see header file for function documentation

  // BAD: unit test fails if I call TGeoGlobalMagField::Instance()
  // at this point. Something goes wrong in the cleaning of the global
  // ROOT onject
  /*
  if (TGeoGlobalMagField::Instance()->GetField()) {
    // everything set, but think about storing the currents for
    // a cross-check
    return 0;
  }
  */

  // The magnetic field is initialized once at the start
  // of run. The fields are supposed to be constant througout the
  // data taking of one run. The run is aborted if the changes
  // exceed a certain limit.
  AliGRPManager grpman;
  if (grpman.ReadGRPEntry() && grpman.SetMagField()) {
    return 0;
  }

  return -ENOENT;
}

AliHLTUInt64_t AliHLTMiscImplementation::GetTriggerMask(AliRawReader* rawReader) const
{
  // see header file for function documentation
  if (!rawReader) return 0;
  AliHLTUInt64_t trgMask=0;
  if (rawReader) {
    const UInt_t* pattern=rawReader->GetTriggerPattern();
    trgMask=pattern[1]&0xfffffff; // 28 upper bits of the 50 bit mask
    trgMask<<=32;
    trgMask|=pattern[0]; // 32 lower bits of the mask
  }
  return trgMask;
}
