// $Id$

//**************************************************************************
//* This file is property of and copyright by the                          * 
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
/// @date   2009-07-07
/// @brief  Implementation of various glue functions implemented in dynamically
///         loaded libraries

#include "AliHLTMiscImplementation.h"
#include "AliHLTLogging.h"
#include "AliHLTErrorGuard.h"
#include "AliLog.h"
#include "AliCDBManager.h"
#include "AliCDBStorage.h"
#include "AliCDBEntry.h"
#include "AliGRPManager.h"
#include "AliRawReader.h"
#include "AliRawEventHeaderBase.h"
#include "AliTracker.h"
#include "AliESDtrack.h"
#include "AliESDHLTDecision.h"
#include "TGeoGlobalMagField.h"
#include "AliHLTGlobalTriggerDecision.h"
#include "TClass.h"
#include "TStreamerInfo.h"
#include "TObjArray.h"
#include <stdexcept>

/** ROOT macro for the implementation of ROOT specific class methods */
ClassImp(AliHLTMiscImplementation);

AliHLTMiscImplementation::AliHLTMiscImplementation()
{
}

AliHLTMiscImplementation::~AliHLTMiscImplementation()
{
  // see header file for function documentation
}

int AliHLTMiscImplementation::InitCDB(const char* cdbpath, const char* cdbsnapshot)
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
  if (cdbsnapshot != NULL) {
    printf("Running in snapshot mode: %s\n", cdbsnapshot);
    pCDB->SetSnapshotMode(cdbsnapshot);
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

AliCDBEntry* AliHLTMiscImplementation::LoadOCDBEntry(const char* path, int runNo) const
{
  // see header file for function documentation
  if (!path) return NULL;
  AliHLTLogging log;

  AliCDBManager* man = AliCDBManager::Instance();
  if (!man) {
    ALIHLTERRORGUARD(1, "failed to access CDB manager");
    return NULL;
  }
  if (runNo<0) runNo=man->GetRun();

  AliCDBEntry* entry=NULL;
  try {
    // exceptions for the loading of OCDB objects have been introduced in r61012 on
    // Feb 20 2013. This allows to reduce this function to try and catch of AliCDBManager::Get
    entry=man->Get(path, runNo);
  }
  catch (...) {
    // ignore the exception, missing object can be a valid condition, and is handled
    // downstream
  }
  return entry;
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
  // check uses AliCDBStorage::GetLatestVersion, which is the only method
  // to check the existence without risking AliFatal

  int iResult=0;
  if (!pMap) return -EINVAL;
  AliHLTLogging log;

  AliCDBManager* man = AliCDBManager::Instance();
  if (!man) {
    ALIHLTERRORGUARD(1, "failed to access CDB manager");
    return -ENOSYS;
  }
  Int_t runNo = GetCDBRunNo();

  TIterator* next=pMap->MakeIterator();
  if (!next) return -EFAULT;

  TObject* pEntry=NULL;
  while ((pEntry=next->Next())) {
    try {
      // exceptions for the loading of OCDB objects have been introduced in r61012 on
      // Feb 20 2013. This allows to reduce this function to try and catch of AliCDBManager::Get
      man->Get(pEntry->GetName(), runNo);
    }
    catch (...) {
      // find out the storage for more details in the error message
      AliCDBStorage* pStorage=NULL;
      const char* uri=man->GetURI(pEntry->GetName());
      if (uri) {
	pStorage = AliCDBManager::Instance()->GetStorage(uri);
      }

      log.Logging(kHLTLogError, "AliHLTMiscImplementation::CheckOCDBEntries", "CDB handling", "can not find required OCDB object %s for run number %d in storage %s", pEntry->GetName(), runNo, (pStorage!=NULL?pStorage->GetURI().Data():"<unavailable>"));
      iResult=-ENOENT;
    }
  }
  delete next;
  next=NULL;

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

AliHLTTriggerMask_t AliHLTMiscImplementation::GetTriggerMask(AliRawReader* rawReader) const
{
  // see header file for function documentation
  if (!rawReader) return 0;
  AliHLTTriggerMask_t trgMask=0;
  if (rawReader) {
    const UInt_t* pattern=rawReader->GetTriggerPattern();
    if(rawReader->GetVersion()==3){
      trgMask=pattern[3];
      trgMask<<=32;
      trgMask|=pattern[2];
      trgMask<<=32;
    }
    trgMask|=pattern[1];
    trgMask<<=32;
    trgMask|=pattern[0]; // 32 lower bits of the mask
  }
  return trgMask;
}

AliHLTUInt32_t AliHLTMiscImplementation::GetTimeStamp(AliRawReader* rawReader) const
{
  // extract time stamp of the event from the event header
  if (!rawReader) return kMaxUInt;
  const AliRawEventHeaderBase* eventHeader = rawReader->GetEventHeader();
  if (!eventHeader) return kMaxUInt;
  return eventHeader->Get("Timestamp");
}

AliHLTUInt32_t AliHLTMiscImplementation::GetEventType(AliRawReader* rawReader) const
{
  // extract event type from the event header
  if (!rawReader) return 0;
  const AliRawEventHeaderBase* eventHeader = rawReader->GetEventHeader();
  if (!eventHeader) return 0;
  UInt_t daqType = eventHeader->Get("Type");
  switch(daqType){
  case AliRawEventHeaderBase::kStartOfRun:
  case AliRawEventHeaderBase::kStartOfData:
    return gkAliEventTypeStartOfRun;

  case AliRawEventHeaderBase::kEndOfRun:
  case AliRawEventHeaderBase::kEndOfData:
    return gkAliEventTypeEndOfRun;

  case AliRawEventHeaderBase::kPhysicsEvent:
    return gkAliEventTypeData;

  case AliRawEventHeaderBase::kCalibrationEvent:
    return gkAliEventTypeCalibration;

  case AliRawEventHeaderBase::kFormatError:
    return gkAliEventTypeCorruptID;

  case AliRawEventHeaderBase::kSystemSoftwareTriggerEvent:
  case AliRawEventHeaderBase::kDetectorSoftwareTriggerEvent:
    return gkAliEventTypeSoftware;

    // TODO: Sync Event Type not implemented!
    //case AliRawEventHeaderBase::kSyncEvent:
  }
  return gkAliEventTypeUnknown;
}

const char* AliHLTMiscImplementation::GetBeamTypeFromGRP() const
{
  // get beam type from GRP
  return NULL;
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
  return AliESDHLTDecision::Class();
}

int AliHLTMiscImplementation::Copy(const AliHLTGlobalTriggerDecision* pDecision, TObject* object) const
{
  // Copy HLT global trigger decision to AliESDHLTDecision container
  if (!pDecision || !object) return -EINVAL;
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
	  TObjArray* pInfos=const_cast<TObjArray*>(pClass->GetStreamerInfos());
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

int AliHLTMiscImplementation::MergeStreamerInfo(TObjArray* tgt, const TObjArray* src, int iVerbosity) const
{
  /// merge streamer info entries from source array to target array
  /// return 1 if target array has been changed

  // add all existing infos if not existing in the current one, or having
  // different class version
  int iResult=0;
  if (!tgt || !src) return -EINVAL;

  {
    // check if all infos from the existing entry are in the new entry and with
    // identical class version
    TIter next(src);
    TObject* nextobj=NULL;
    while ((nextobj=next())) {
      TStreamerInfo* srcInfo=dynamic_cast<TStreamerInfo*>(nextobj);
      if (!srcInfo) continue;
      TString srcInfoName=srcInfo->GetName();

      int i=0;
      for (; i<tgt->GetEntriesFast(); i++) {
	if (tgt->At(i)==NULL) continue;
	if (srcInfoName.CompareTo(tgt->At(i)->GetName())!=0) continue;
	// TODO: 2010-08-23 some more detailed investigation is needed.
	// Structures used for data exchange, e.g. AliHLTComponentDataType
	// or AliHLTEventDDLV1 do not have a class version, but need to be stored in the
	// streamer info. Strictly speaking not, because those structures are not supposed
	// to be changed at all, so they should be the same in all versions in the future.
	// There has been a problem with detecting whether the streamer info is already in
	// the target array if the srcInfo has class version -1. As it just concerns
	// structures not going to be changed we can safely skip checking the class version,
	// as long as the entry is already in the target streamer infos it does not need
	// to be copied again.
	if (srcInfo->GetClassVersion()<0) break;
	TStreamerInfo* tgtInfo=dynamic_cast<TStreamerInfo*>(tgt->At(i));
	if (tgtInfo && tgtInfo->GetClassVersion()==srcInfo->GetClassVersion()) break;
      }
      if (i<tgt->GetEntriesFast()) continue;

      iResult=1;
      if (iVerbosity>0) {
	AliInfo(Form("adding streamer info for class %s version %d", srcInfoName.Data(), srcInfo->GetClassVersion()));
      }
      tgt->Add(srcInfo);
    }
  }

  return iResult;
}

void AliHLTMiscImplementation::SetAliESDtrackOnlineModeFlag(bool mode) const
{
  /// set the online mode flag of AliESDtrack
  AliESDtrack::OnlineMode(mode);
}

bool AliHLTMiscImplementation::GetAliESDtrackOnlineModeFlag() const
{
  /// get status of the online mode flag of AliESDtrack
  return AliESDtrack::OnlineMode();
}
