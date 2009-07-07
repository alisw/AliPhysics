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
#include "AliCDBEntry.h"

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
      pCDB->SetDefaultStorage(cdbpath);
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
    log.Logging(kHLTLogError, "InitCDB", "CDB handling", "Could not get CDB instance");
  } else {
    pCDB->SetRun(runNo);
  }
  return iResult;
}

AliCDBEntry* AliHLTMiscImplementation::LoadOCDBEntry(const char* path, int runNo, int version, int subVersion)
{
  // see header file for function documentation
  AliCDBManager::Instance()->UnloadFromCache(path);
  return AliCDBManager::Instance()->Get(path, runNo, version, subVersion);
}

TObject* AliHLTMiscImplementation::ExtractObject(AliCDBEntry* entry)
{
  // see header file for function documentation
  if (!entry) return NULL;
  return entry->GetObject();
}

int AliHLTMiscInitCDB(const char* cdbpath)
{
  int iResult=0;
  AliCDBManager* pCDB = AliCDBManager::Instance();
  AliHLTLogging log;
  if (!pCDB) {
    log.Logging(kHLTLogError, "InitCDB", "CDB handling", "Could not get CDB instance");
  } else {
    if (cdbpath && cdbpath[0]!=0) {
      pCDB->SetDefaultStorage(cdbpath);
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

int AliHLTMiscSetCDBRunNo(int runNo)
{
  int iResult=0;
  AliCDBManager* pCDB = AliCDBManager::Instance();
  AliHLTLogging log;
  if (!pCDB) {
    log.Logging(kHLTLogError, "InitCDB", "CDB handling", "Could not get CDB instance");
  } else {
    pCDB->SetRun(runNo);
  }
  return iResult;
}
