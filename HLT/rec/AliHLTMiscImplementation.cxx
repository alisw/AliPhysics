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
#include "AliLog.h"
#include "AliCDBManager.h"
#include "AliCDBStorage.h"
#include "AliCDBEntry.h"
#include "AliGRPManager.h"
#include "AliRawReader.h"
#include "AliTracker.h"
#ifndef HAVE_NOT_ALIESDHLTDECISION
#include "AliESDHLTDecision.h"
#endif //HAVE_NOT_ALIESDHLTDECISION
#include "TGeoGlobalMagField.h"
#include "AliHLTGlobalTriggerDecision.h"
#include "TClass.h"
#include "TStreamerInfo.h"
#include "TObjArray.h"

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

int AliHLTMiscImplementation::GetCDBRunNo() const
{
  // see header file for function documentation
  AliCDBManager* pCDB = AliCDBManager::Instance();
  AliHLTLogging log;
  if (!pCDB) {
    log.Logging(kHLTLogError, "SetCDBRunNo", "CDB handling", "Could not get CDB instance");
  } else {
    return pCDB->GetRun();
  }
  return -1;
}

AliCDBEntry* AliHLTMiscImplementation::LoadOCDBEntry(const char* path, int runNo, int version, int subVersion) const
{
  // see header file for function documentation
  AliCDBStorage* store = AliCDBManager::Instance()->GetDefaultStorage();
  if (!store) {
    return NULL;
  }
  if (version<0) version = store->GetLatestVersion(path, runNo);
  if (subVersion<0) subVersion = store->GetLatestSubVersion(path, runNo, version);
  return store->Get(path, runNo, version, subVersion);
}

TObject* AliHLTMiscImplementation::ExtractObject(AliCDBEntry* entry) const
{
  // see header file for function documentation
  if (!entry) return NULL;
  return entry->GetObject();
}

int AliHLTMiscImplementation::CheckOCDBEntries(const TMap* const pMap) const
{
  // check the availability of the OCDB entry descriptions in the TMap
  //  key : complete OCDB path of the entry
  //  value : auxiliary object - short description
  int iResult=0;
  if (!pMap) return -EINVAL;

  const TMap* pStorages=AliCDBManager::Instance()->GetStorageMap();
  Int_t runNo = GetCDBRunNo();

  TIterator* next=pMap->MakeIterator();
  TObject* pEntry=NULL;
  while ((pEntry=next->Next())) {
    // check if the entry has specific storage
    AliCDBStorage* pStorage=NULL;
    TObject* pStorageId=pStorages->GetValue(pEntry->GetName());
    if (pStorageId) {
      pStorage=AliCDBManager::Instance()->GetStorage(pStorageId->GetName());
    } else {
      pStorage=AliCDBManager::Instance()->GetDefaultStorage();
    }

    if (pStorage->GetLatestVersion(pEntry->GetName(), runNo)<0) {
      AliHLTLogging log;
      log.Logging(kHLTLogError, "CheckOCDBEntries", "CDB handling", "can not find required OCDB object %s for run number %d in storage %s", pEntry->GetName(), runNo, pStorage->GetURI().Data());
      iResult=-ENOENT;
    } else {
      AliHLTLogging log;
      log.Logging(kHLTLogDebug, "CheckOCDBEntries", "CDB handling", "found required OCDB object %s for run number %d in storage %s", pEntry->GetName(), runNo, pStorage->GetURI().Data());
    }
  }

  return iResult;
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

Double_t AliHLTMiscImplementation::GetBz()
{
  // Returns Bz.
  return AliTracker::GetBz();
}

Double_t AliHLTMiscImplementation::GetBz(const Double_t *r)
{
  // Returns Bz (kG) at the point "r" .
  return AliTracker::GetBz(r);
}

void AliHLTMiscImplementation::GetBxByBz(const Double_t r[3], Double_t b[3])
{
  // Returns Bx, By and Bz (kG) at the point "r" .
  return AliTracker::GetBxByBz(r, b);
}

const TClass* AliHLTMiscImplementation::IsAliESDHLTDecision() const
{
  // Return the IsA of the AliESDHLTDecision class
#ifndef HAVE_NOT_ALIESDHLTDECISION
  return AliESDHLTDecision::Class();
#else // HAVE_NOT_ALIESDHLTDECISION
  return NULL;
#endif // HAVE_NOT_ALIESDHLTDECISION
}

int AliHLTMiscImplementation::Copy(const AliHLTGlobalTriggerDecision* pDecision, TObject* object) const
{
  // Copy HLT global trigger decision to AliESDHLTDecision container
  if (!pDecision || !object) return -EINVAL;
#ifndef HAVE_NOT_ALIESDHLTDECISION
  AliESDHLTDecision* pESDHLTDecision=NULL;
  if (object->IsA()==NULL ||
      object->IsA() != AliESDHLTDecision::Class() ||
      (pESDHLTDecision=dynamic_cast<AliESDHLTDecision*>(object))==NULL) {
//     HLTError("can not copy HLT global decision to object of class \"%s\"", 
// 	     object->IsA()?object->IsA()->GetName():"NULL");
    return -EINVAL;
  }

  pESDHLTDecision->~AliESDHLTDecision();
  new (pESDHLTDecision) AliESDHLTDecision(pDecision->Result(), pDecision->GetTitle());

#endif // HAVE_NOT_ALIESDHLTDECISION
  return 0;
}

int AliHLTMiscImplementation::InitStreamerInfos(const char* ocdbEntry) const
{
  // init streamer infos for HLT reconstruction
  // Root schema evolution is not enabled for AliHLTMessage and all streamed objects.
  // Objects in the raw data payload rely on the availability of the correct stream info.
  // The relevant streamer info for a specific run is stored in the OCDB.
  int iResult=0;

  AliCDBEntry* pEntry=LoadOCDBEntry(ocdbEntry);
  TObject* pObject=NULL;
  if (pEntry && (pObject=ExtractObject(pEntry))!=NULL)
    {
    TObjArray* pSchemas=dynamic_cast<TObjArray*>(pObject);
    if (pSchemas) {
      iResult=InitStreamerInfos(pSchemas);
    } else {
      AliError(Form("internal mismatch in OCDB entry %s: wrong class type", ocdbEntry));
    }
  } else {
    AliWarning(Form("missing HLT reco data (%s), skipping initialization of streamer info for TObjects in HLT raw data payload", ocdbEntry));
  }
  return iResult;
}

int AliHLTMiscImplementation::InitStreamerInfos(TObjArray* pSchemas) const
{
  // init streamer infos for HLT reconstruction from an array of TStreamerInfo objects

  for (int i=0; i<pSchemas->GetEntriesFast(); i++) {
    if (pSchemas->At(i)) {
      TStreamerInfo* pSchema=dynamic_cast<TStreamerInfo*>(pSchemas->At(i));
      if (pSchema) {
	int version=pSchema->GetClassVersion();
	TClass* pClass=TClass::GetClass(pSchema->GetName());
	if (pClass) {
	  if (pClass->GetClassVersion()==version) {
	    AliDebug(0,Form("skipping schema definition %d version %d to class %s as this is the native version", i, version, pSchema->GetName()));
	    continue;
	  }
	  TObjArray* pInfos=pClass->GetStreamerInfos();
	  if (pInfos /*&& version<pInfos->GetEntriesFast()*/) {
	    if (pInfos->At(version)==NULL) {
	      TStreamerInfo* pClone=(TStreamerInfo*)pSchema->Clone();
	      if (pClone) {
		pClone->SetClass(pClass);
		pClone->BuildOld();
		pInfos->AddAtAndExpand(pClone, version);
		AliDebug(0,Form("adding schema definition %d version %d to class %s", i, version, pSchema->GetName()));
	      } else {
		AliError(Form("skipping schema definition %d (%s), unable to create clone object", i, pSchema->GetName()));
	      }
	    } else {
	      TStreamerInfo* pInfo=dynamic_cast<TStreamerInfo*>(pInfos->At(version));
	      if (pInfo && pInfo->GetClassVersion()==version) {
		AliDebug(0,Form("schema definition %d version %d already available in class %s, skipping ...", i, version, pSchema->GetName()));
	      } else {
		AliError(Form("can not verify version for already existing schema definition %d (%s) version %d: version of existing definition is %d", i, pSchema->GetName(), version, pInfo?pInfo->GetClassVersion():-1));
	      }
	    }
	  } else {
	    AliError(Form("skipping schema definition %d (%s), unable to set version %d in info array of size %d", i, pSchema->GetName(), version, pInfos?pInfos->GetEntriesFast():-1));
	  }
	} else {
	  AliError(Form("skipping schema definition %d (%s), unable to find class", i, pSchema->GetName()));
	}
      } else {
	AliError(Form("skipping schema definition %d, not of TStreamerInfo", i));
      }
    }
  }

  return 0;
}
