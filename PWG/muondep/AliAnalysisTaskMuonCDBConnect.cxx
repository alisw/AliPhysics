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

// ROOT includes
#include <TGeoManager.h>
#include <TGeoGlobalMagField.h>

// STEER includes
#include "AliCDBManager.h"
#include "AliCDBStorage.h"
#include "AliGeomManager.h"
#include "AliLog.h"

// ANALYSIS includes
#include "AliAnalysisTaskMuonCDBConnect.h"

// MUON includes
#include "AliMUONCDB.h"
#include "AliMUONESDInterface.h"
#include "AliMUONVTrackReconstructor.h"
#include "AliMUONRecoParam.h"

// MUON mapping includes
#include "AliMpSegmentation.h"
#include "AliMpDDLStore.h"
#include "AliMpManuStore.h"

ClassImp(AliAnalysisTaskMuonCDBConnect)

//________________________________________________________________________
AliAnalysisTaskMuonCDBConnect::AliAnalysisTaskMuonCDBConnect() :
AliAnalysisTaskSE(),
fDefaultStorage(""),
fAlignStorage(""),
fAlignVersion(-1),
fAlignSubVersion(-1),
fRecoParamStorage(""),
fLoadMagField(kFALSE),
fLoadGeometry(kFALSE),
fLoadMapping(kFALSE),
fSegOnly(kFALSE),
fOCDBSet(kFALSE)
{
  /// Default constructor
}

//________________________________________________________________________
AliAnalysisTaskMuonCDBConnect::AliAnalysisTaskMuonCDBConnect(const char *name) :
AliAnalysisTaskSE(name),
fDefaultStorage("raw://"),
fAlignStorage(""),
fAlignVersion(-1),
fAlignSubVersion(-1),
fRecoParamStorage(""),
fLoadMagField(kFALSE),
fLoadGeometry(kFALSE),
fLoadMapping(kFALSE),
fSegOnly(kFALSE),
fOCDBSet(kFALSE)
{
  /// Constructor
}

//________________________________________________________________________
Bool_t AliAnalysisTaskMuonCDBConnect::UserNotify()
{
  /// setup OCDB default storage
  
  // do it only once
  if (fOCDBSet) return kTRUE;
  
  // set OCDB location
  AliCDBManager* cdbm = AliCDBManager::Instance();
  if (!cdbm->IsDefaultStorageSet()) cdbm->SetDefaultStorage(fDefaultStorage.Data());
  else printf("MuonCDBConnect: CDB default storage already set. Do not change it.\n");
  
  return kTRUE;
  
}

//________________________________________________________________________
void AliAnalysisTaskMuonCDBConnect::NotifyRun()
{
  /// setup OCDB path and load requested data
  
  // do it only once
  if (fOCDBSet) return;
  
  // set run number
  AliCDBManager* cdbm = AliCDBManager::Instance();
  if (cdbm->GetRun() < 0) cdbm->SetRun(fCurrentRunNumber);
  else printf("MuonCDBConnect: run number already set. Do not change it.\n");
  
  // set specific storage for MUON alignment
  if (!fAlignStorage.IsNull() || fAlignVersion >= 0 || fAlignSubVersion >= 0) {
    if (!AliGeomManager::GetGeometry()) {
      if (fAlignStorage != "none") {
        if (fAlignStorage.IsNull()) cdbm->SetSpecificStorage("MUON/Align/Data", cdbm->GetDefaultStorage()->GetURI().Data(), fAlignVersion, fAlignSubVersion);
        else cdbm->SetSpecificStorage("MUON/Align/Data", fAlignStorage.Data(), fAlignVersion, fAlignSubVersion);
      }
    } else printf("MuonCDBConnect: geometry already loaded. Do not change MUON align storage.\n");
  }
  
  // set specific storage for MUON recoParam
  if (!fRecoParamStorage.IsNull()) {
    if (!AliMUONESDInterface::GetTracker() || !AliMUONESDInterface::GetTracker()->GetRecoParam())
      cdbm->SetSpecificStorage("MUON/Calib/RecoParam",fRecoParamStorage.Data());
    else printf("MuonCDBConnect: MUON recoParam already loaded. Do not change MUON recoParam storage.\n");
  }
  
  // load magnetic field
  if (fLoadMagField) {
    if (TGeoGlobalMagField::Instance()->GetField()) printf("MuonCDBConnect: magnetic field already loaded. Do not change it.\n");
    else if (!AliMUONCDB::LoadField()) AliFatal("loading of magnetic field failed");
  }
  
  // load geometry
  if (fLoadGeometry) {
    if (AliGeomManager::GetGeometry()) printf("MuonCDBConnect: geometry already loaded. Do not change it.\n");
    else {
      AliGeomManager::LoadGeometry();
      if (!AliGeomManager::GetGeometry() || (fAlignStorage != "none" && !AliGeomManager::ApplyAlignObjsFromCDB("MUON")))
        AliFatal("loading of geometry failed");
    }
  }
  
  // load mapping
  if (fLoadMapping) {
    if ((!fSegOnly && AliMpDDLStore::Instance(kFALSE) && AliMpManuStore::Instance(kFALSE)) ||
        (fSegOnly && AliMpSegmentation::Instance(kFALSE))) printf("MuonCDBConnect: mapping already loaded. Do not change it.\n");
    else if (!AliMUONCDB::LoadMapping(fSegOnly)) AliFatal("loading of mapping failed");
  }
  
  fOCDBSet = kTRUE;
  
}

