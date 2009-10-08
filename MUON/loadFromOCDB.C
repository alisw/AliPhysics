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

// $Id$

/// \ingroup macros
/// \file loadFromOCDB.C
/// \brief Load magnetic field, mapping and reconstruction parameters from OCDB
/// \author Philippe Pillot, SUBATECH

#include <TString.h>
#include <TObjArray.h>
#include <TGeoGlobalMagField.h>
#include <TGrid.h>

#include "AliCDBManager.h"
#include "AliCDBEntry.h"
#include "AliGRPManager.h"

#include "AliMpCDB.h"

#include "AliMUONRecoParam.h"

// prototypes
Bool_t setupOCDB(const TString ocdbPath = "local://$ALICE_ROOT/OCDB", Int_t runNumber = 0);
Bool_t checkOCDB();
Bool_t loadField();
Bool_t loadMapping();
AliMUONRecoParam* loadRecoParam();

//----------------------------------------------------------------------------------------
AliMUONRecoParam* loadFromOCDB(const TString ocdbPath = "local://$ALICE_ROOT/OCDB", Int_t runNumber = 0,
			       Bool_t ldField = kTRUE, Bool_t ldMapping = kTRUE, Bool_t ldRecoParam = kTRUE)
{
  /// Set the default OCDB path and run number (overload the existing ones if already set).
  /// Load magnetic field, mapping and reconstruction parameters from OCDB as requested by the flags.
  /// Return the AliMUONRecoParam object if loaded (0x0 otherwise).
  /// Exit in case of any failure.
  
  AliMUONRecoParam* recoParam = 0x0;
  
  if (!setupOCDB(ocdbPath, runNumber)) exit(-1);
  
  if (ldField && !loadField()) exit(-1);
  
  if (ldMapping && !loadMapping()) exit(-1);
  
  if (ldRecoParam && !(recoParam = loadRecoParam())) exit(-1);
  
  return recoParam;
  
}

//----------------------------------------------------------------------------------------
Bool_t setupOCDB(const TString ocdbPath, Int_t runNumber)
{
  /// Set the default OCDB storage and run number (overload the existing ones if already set).
  
  Info("setupOCDB","setting OCDB path...");
  
  if (runNumber < 0) {
    Error("setupOCDB", "runNumber must be a positive value");
    return kFALSE;
  }
  
  if (ocdbPath.BeginsWith("alien://")) TGrid::Connect("alien://");
  
  AliCDBManager* man = AliCDBManager::Instance();
  
  man->SetDefaultStorage(ocdbPath.Data());
  if (!man->IsDefaultStorageSet()) return kFALSE;
  
  man->SetRun(runNumber);
  
  return kTRUE;
  
}

//----------------------------------------------------------------------------------------
Bool_t checkOCDB()
{
  /// Check that OCDB path and run number are properly set
  
  AliCDBManager* man = AliCDBManager::Instance();
  if (!man->IsDefaultStorageSet() || man->GetRun() < 0) {
    Error("checkOCDB", "OCDB path and runNumber must be properly set");
    return kFALSE;
  }
  
  return kTRUE;
  
}

//----------------------------------------------------------------------------------------
Bool_t loadField()
{
  /// Load magnetic field (existing field will be deleted).
  /// OCDB path is supposed to be set.
  
  Info("loadField","Loading field map from GRP...");
  
  if (!checkOCDB()) return kFALSE;
  
  AliGRPManager grpMan;
  
  // in case it has already been set
  TGeoGlobalMagField::Instance()->Unlock();
  
  if (!grpMan.ReadGRPEntry() || !grpMan.SetMagField()) {
    Error("loadField", "failed to load magnetic field from OCDB");
    return kFALSE;
  }
  
  return kTRUE;
  
}

//----------------------------------------------------------------------------------------
Bool_t loadMapping()
{
  /// Load mapping (existing mapping will be unloaded).
  /// OCDB path is supposed to be set.
  
  Info("loadMapping","Loading mapping from OCDB...");
  
  if (!checkOCDB()) return kFALSE;
  
  // in case it has already been set
  AliMpCDB::UnloadAll();
  
  if (!AliMpCDB::LoadAll(kTRUE)) {
    Error("loadMapping","failed to load mapping from OCDB");
    return kFALSE;
  }
  
  return kTRUE;
  
}

//----------------------------------------------------------------------------------------
AliMUONRecoParam* loadRecoParam()
{
  /// Load and return reconstruction parameters.
  /// OCDB path is supposed to be set.
  
  Info("loadRecoParam","Loading RecoParam from OCDB...");
  
  if (!checkOCDB()) return kFALSE;
  
  AliMUONRecoParam* recoParam = 0x0;
  AliCDBEntry* entry = AliCDBManager::Instance()->Get("MUON/Calib/RecoParam");
  
  if(entry) {
    
    // load recoParam according OCDB content (single or array)
    if (!(recoParam = dynamic_cast<AliMUONRecoParam*>(entry->GetObject()))) {
      
      TObjArray* recoParamArray = static_cast<TObjArray*>(entry->GetObject());
      
      for(Int_t i = 0; i < recoParamArray->GetEntriesFast(); i++) {
	recoParam = static_cast<AliMUONRecoParam*>(recoParamArray->UncheckedAt(i));
	if (recoParam->IsDefault()) break;
	recoParam = 0x0;
      }
      
    }
    
  }
  
  if (!recoParam) Error("loadRecoParam", "failed to load RecoParam from OCDB");
  
  return recoParam;
  
}

